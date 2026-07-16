using DataBrowser
using DataBrowserAnnotations
using Test

const COLLECTION_ID_INDEX = DataBrowserAPI.ItemIndex
const COLLECTION_ID_BROWSER = DataBrowserGUI.Browser

struct DurableCollectionSource <: AbstractDataSource
    root_path::String
    label_suffix::String
end

struct DurableCollectionSourceItem <: AbstractDataSourceItem
    id::String
end

struct DurableCollection <: AbstractCollection
    key::Int
    shown::String
end

struct DurableCollectionItem <: AbstractDataItem
    id::String
    path::Vector{AbstractCollection}
end

DataBrowser.source_id(source::DurableCollectionSource) = source.root_path
DataBrowser.source_label(::DurableCollectionSource) = "durable collection source"
DataBrowser.source_items(::DurableCollectionSource; kwargs...) =
    [DurableCollectionSourceItem("source")]
DataBrowser.source_item_id(item::DurableCollectionSourceItem) = item.id
DataBrowser.source_item_label(item::DurableCollectionSourceItem) = item.id
DataBrowser.fingerprint(::DurableCollectionSourceItem) = "stable source"

DataBrowser.id(collection::DurableCollection) = collection.key
DataBrowser.label(collection::DurableCollection) = collection.shown

DataBrowser.id(item::DurableCollectionItem) = item.id
DataBrowser.collection(item::DurableCollectionItem) = item.path

function DataBrowser.data_items(
    ::Project,
    source::DurableCollectionSource,
    ::DurableCollectionSourceItem,
)
    parent = DurableCollection(1, "parent $(source.label_suffix)")
    return [
        DurableCollectionItem(
            "durable-1",
            AbstractCollection[parent, DurableCollection(11, "same label")],
        ),
        DurableCollectionItem(
            "durable-2",
            AbstractCollection[parent, DurableCollection(12, "same label")],
        ),
    ]
end

@testset "collection IDs survive a clean cache rebuild" begin
    root_path = mktempdir()
    project = define_project("CollectionIdPersistence")
    first = open_workspace(
        project,
        DurableCollectionSource(root_path, "before");
        cache=true,
        background_processing=true,
    )

    cache_path = first.cache.identity.cache_path
    annotation_root = dirname(cache_path)
    durable_id = ""
    annotation_key = ""
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(first)
        collections = first.index.collections
        first_record = first.index.items["durable-1"]
        second_record = first.index.items["durable-2"]
        first_collection_record = collections.records[first_record.collection_key]
        second_collection_record = collections.records[second_record.collection_key]

        @test first_collection_record.label == second_collection_record.label == "same label"
        @test first_collection_record.id != second_collection_record.id
        durable_id = first_collection_record.id
        annotation_key = COLLECTION_ID_BROWSER._collection_annotation_key(
            collections, first_collection_record)
        tag_state = DataBrowserAnnotations.Tags.TagState(
            [DataBrowserAnnotations.Tags.TagDef("review", (0x30, 0xc0, 0xff), 50)],
            Dict(annotation_key => Set(["review"])),
        )
        DataBrowserAnnotations.Tags.save(annotation_root, tag_state)
        DataBrowserAnnotations.Notes.write_section!(
            annotation_root, annotation_key, "persists across rebuild")

        view = COLLECTION_ID_BROWSER.PersistedProjectView(
            project=project.name,
            tree=COLLECTION_ID_BROWSER.PersistedTreeView(
                expanded=[durable_id],
                selected=[durable_id],
            ),
            items=COLLECTION_ID_BROWSER.PersistedItemsView(
                selected=[first_record.id],
            ),
        )
        COLLECTION_ID_BROWSER._save_project_view(root_path, view)
    finally
        close_workspace!(first)
    end

    rm(cache_path; force=true)
    rm(cache_path * ".wal"; force=true)

    reopened = open_workspace(
        project,
        DurableCollectionSource(root_path, "after");
        rebuild=true,
        cache=true,
        background_processing=true,
    )
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(reopened)
        collections = reopened.index.collections
        first_record = reopened.index.items["durable-1"]
        second_record = reopened.index.items["durable-2"]
        first_collection_record = collections.records[first_record.collection_key]
        second_collection_record = collections.records[second_record.collection_key]

        @test first_collection_record.id == durable_id
        root_key = only(COLLECTION_ID_INDEX.sorted_child_keys(collections, nothing))
        @test collections.records[root_key].label == "parent after"

        state = COLLECTION_ID_BROWSER.BrowserState(workspace=reopened)
        COLLECTION_ID_BROWSER._load_tag_state_for_root!(state, annotation_root)
        view = COLLECTION_ID_BROWSER._load_project_view(root_path)
        COLLECTION_ID_BROWSER._apply_project_view!(state, view)
        selected_collections, selected_items, _ =
            COLLECTION_ID_BROWSER._project_visible_selection(state)

        @test state.expanded_collection_ids == [durable_id]
        @test only(selected_collections).id == durable_id
        @test only(selected_items).id == "durable-1"
        @test state.tag_state.assignments[annotation_key] == Set(["review"])
        @test !haskey(
            state.tag_state.assignments,
            COLLECTION_ID_BROWSER._collection_annotation_key(collections, second_collection_record),
        )
        @test DataBrowserAnnotations.Notes.read_section(annotation_root, annotation_key) ==
            "persists across rebuild"
    finally
        close_workspace!(reopened)
        rm(root_path; force=true, recursive=true)
        rm(annotation_root; force=true, recursive=true)
    end
end
