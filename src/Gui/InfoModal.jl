import CImGui as ig
import CImGui.CSyntax: @c

using ..Projects:
    kind_label
import ..Workspace

"""Render collection and item details for the visible workspace selection."""
function render_info_window(state::BrowserState)::Nothing
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        if ig.Begin("Information Panel")
            ig.TextDisabled("Open a project folder to inspect items")
        end
        ig.End()
        return nothing
    end
    selected_collections, selected_items, selected_path =
        _project_visible_selection(state)
    if ig.Begin("Information Panel")
        flags = ig.ImGuiTableFlags_Borders | ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_ScrollY
        ig.BeginTable("info_cols", 2, flags)
        ig.TableSetupColumn("Collection")
        ig.TableSetupColumn("Item")
        ig.TableHeadersRow()
        ig.TableNextRow()
        ig.TableNextColumn()

        if length(selected_collections) == 1
            collection_node = selected_collections[1]
            sel_name = join(selected_path, "/")
            ig.Text("Location: $sel_name")
            ig.Separator()
            collection_metadata = merge(collection_node.metadata, collection_node.analysis)
            if !isempty(collection_metadata)
                ig.Text("Collection metadata")
                for k in sort!(collect(keys(collection_metadata)); by=String)
                    ig.BulletText("$(k): $(collection_metadata[k])")
                end
            else
                ig.TextDisabled("No collection metadata found")
            end
        elseif isempty(selected_collections)
            ig.TextDisabled("Select a collection to see details")
        else
            ig.TextDisabled("Select a single collection to see details")
        end

        ig.TableNextColumn()
        if length(selected_items) == 1
            m = selected_items[1]
            # The delivered metadata dict (inherited ⊕ entries ⊕ computed layers), exactly what
            # project callbacks and views receive.
            delivered = Workspace.delivered_metadata(workspace, m)
            ig.Text("Title: $(m.item_label)")
            ig.Separator()
            kind_text = kind_label(workspace.project, m.kind)
            ig.BulletText("Type: $(kind_text)")
            ig.BulletText("Timestamp: $(m.source_item_timestamp)")
            if m.source_item_path !== nothing
                ig.BulletText("Source item:")
                ig.SameLine()
                ig.TextLinkOpenURL(m.source_item_id, m.source_item_path)
            else
                ig.BulletText("Source item: $(m.source_item_id)")
            end
            ig.Separator()
            if !isempty(delivered)
                ig.Text("Metadata")
                for k in sort!(collect(keys(delivered)); by=String)
                    ig.BulletText("$(k) = $(delivered[k])")
                end
            else
                ig.TextDisabled("No metadata")
            end

        elseif isempty(selected_items)
            ig.TextDisabled("Select an item to view details")
        else
            ig.TextDisabled("Select a single item to view details")
        end
        ig.EndTable()
    end
    ig.End()
    return nothing
end

"""Tell the user when the current source did not provide collection parameters."""
function render_collection_metadata_modal(state::BrowserState)::Nothing
    # Dear ImGui cannot hold two sibling root modals open: their per-frame OpenPopup calls replace
    # each other in the popup stack, leaving an invisible modal that swallows all input. Yield to
    # the cache-rebuild modal and appear after it is dismissed.
    state.cache_rebuild_modal && return nothing
    workspace = state.workspace
    workspace isa Workspace.Workspace || return nothing
    # Reset dismissal when root path changes
    current_root = workspace.cache.identity.source_id
    if state.modal_root_path != current_root
        state.modal_root_path = current_root
        state.collection_metadata_modal = true
    end
    # always center
    center = ig.ImVec2(0.5, 0.5)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))

    if state.collection_metadata_modal &&
       !workspace.index.hierarchy.has_collection_metadata
        ig.OpenPopup("Collection Parameters Missing")
    end

    opened = state.collection_metadata_modal

    if @c ig.BeginPopupModal("Collection Parameters Missing", &opened, ig.ImGuiWindowFlags_AlwaysAutoResize)
        ig.Text("No collection parameters were provided by this source.")
        ig.Separator()
        ig.TextWrapped(
            "The browser will still work, but collection-level fields such as area, thickness, " *
            "or layout coordinates will be empty until the source supplies them.",
        )
        ig.Spacing()
        if ig.Button("Got it")
            opened = false
            ig.CloseCurrentPopup()
        end
        ig.EndPopup()
    end
    state.collection_metadata_modal = opened
    return nothing
end

"""Ask before discarding an old generated cache whose schema cannot be opened."""
function render_cache_rebuild_modal(state::BrowserState)::Nothing
    state.cache_rebuild_modal || return nothing

    center = ig.ImVec2(0.5, 0.5)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))
    ig.OpenPopup("Cache Rebuild Required")

    opened = state.cache_rebuild_modal
    if @c ig.BeginPopupModal("Cache Rebuild Required", &opened, ig.ImGuiWindowFlags_AlwaysAutoResize)
        ig.TextWrapped("This project cache was made by an older MeasurementBrowser schema.")
        ig.TextWrapped(
            "The workspace is already open without disk persistence. Source data is safe; only " *
            "the generated cache file would be discarded.",
        )
        ig.Separator()
        ig.TextWrapped(state.cache_rebuild_error)
        ig.Spacing()

        if ig.Button("Continue without disk cache")
            opened = false
            state.cache_rebuild_path = ""
            state.cache_rebuild_project = nothing
            state.cache_rebuild_error = ""
            ig.CloseCurrentPopup()
        end
        ig.SameLine()
        if ig.Button("Discard cache and rebuild")
            path = state.cache_rebuild_path
            project = state.cache_rebuild_project
            persist = state.cache_rebuild_persist
            project === nothing && error("Cannot rebuild cache because no project was selected")
            opened = false
            state.cache_rebuild_modal = false
            state.cache_rebuild_path = ""
            state.cache_rebuild_project = nothing
            state.cache_rebuild_error = ""
            ig.CloseCurrentPopup()
            _open_project_path!(state, path; project, persist, rebuild_cache=true)
        end
        ig.EndPopup()
    end
    state.cache_rebuild_modal = opened
    return nothing
end
