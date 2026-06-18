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
            item_vec = selected_collections[1].items
            sel_name = join(selected_path, "/")
            ig.Text("Location: $sel_name")
            ig.Separator()
            if !isempty(item_vec)
                collection_meta = first(item_vec).collection_metadata
                if !isempty(collection_meta)
                    ig.Text("Collection metadata")
                    for (k, v) in collection_meta
                        ig.BulletText("$(k): $(v)")
                    end
                else
                    ig.TextDisabled("No metadata parameters found")
                end
            end
        elseif isempty(selected_collections)
            ig.TextDisabled("Select a collection to see details")
        else
            ig.TextDisabled("Select a single collection to see details")
        end

        ig.TableNextColumn()
        if length(selected_items) == 1
            m = selected_items[1]
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
            if !isempty(m.parameters)
                ig.Text("Parameters")
                for k in sort!(collect(keys(m.parameters)); by=String)
                    v = m.parameters[k]
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No parameters extracted")
            end

            if !isempty(m.stats)
                ig.Separator()
                ig.Text("Statistics")
                for k in sort!(collect(keys(m.stats)); by=String)
                    v = m.stats[k]
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No stats computed")
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

"""Tell the user when the current source did not provide collection metadata."""
function render_collection_metadata_modal(state::BrowserState)::Nothing
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

    # Show modal if: missing metadata and user hasn't dismissed it this scan
    if state.collection_metadata_modal &&
       !workspace.index.hierarchy.has_collection_metadata
        ig.OpenPopup("Collection Metadata Missing")
    end

    opened = state.collection_metadata_modal

    if @c ig.BeginPopupModal("Collection Metadata Missing", &opened, ig.ImGuiWindowFlags_AlwaysAutoResize)
        ig.Text("No collection metadata was provided by this source.")
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
