module Browser

using Annotations
using DataFrames: DataFrame
using Printf
using Statistics: mean

import GLFW
using GLMakie
import GLMakie.Makie as Makie
import CImGui as ig
import CImGui.CSyntax: @c
import ModernGL as gl
using NativeFileDialog: pick_folder, save_file
using TOML

using ..Project:
    AbstractProject,
    DEFAULT_PROJECT,
    PROJECTS,
    display_label,
    kind_label,
    load_source_data,
    project_description,
    project_name
using ..MeasurementIndex:
    HierarchyNode,
    MeasurementInfo,
    SourceScan,
    children,
    device_path_key,
    device_path_tuple,
    index_source_file,
    measurement_timestamp_key
using ..Cache:
    ProjectCacheIdentity,
    ProjectCacheStatus
import ..Workspace
using ..Workspace:
    WorkspaceProgress,
    cache_work_running,
    cancel_cache!,
    cancel_scan!,
    close_workspace!,
    open_workspace,
    poll_workspace!,
    scan_source!,
    source_scan_running,
    track_task!,
    update_cache!
using ..Visualization:
    PlotKind,
    debug_plot,
    plot_data!,
    plot_kinds,
    setup_plot
using ..MeasurementBrowser:
    FigureScriptExistsError,
    FigureScriptIOError,
    FigureScriptInferenceProfile,
    FigureScriptProfiledError,
    FigureScriptValidationError,
    NamedMeasurementGroup,
    _FigureScriptFactIndex,
    _build_figure_script_fact_index,
    _infer_measurement_group_profiled,
    _matching_measurements,
    _validate_named_measurement_groups,
    figure_script_path,
    write_figure_script

include("Browser/MakieIntegration.jl")
using .MakieImguiIntegration: MakieFigure

include("Gui/BrowserState.jl")
include("Gui/BadAndStyling.jl")
include("Gui/BrowserOperations.jl")
include("Gui/PlotPanel.jl")
include("Gui/Persistence.jl")
include("Gui/TreePanel.jl")
include("Gui/InfoModal.jl")
include("Gui/Layout.jl")

export start_browser

end
