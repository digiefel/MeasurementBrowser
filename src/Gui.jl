import GLFW
using GLMakie
import GLMakie.Makie as Makie
import CImGui as ig
import CImGui.CSyntax: @c
import ModernGL as gl
using Printf

using Statistics: mean

include("MakieIntegration.jl")
using .MakieImguiIntegration

using TOML
using NativeFileDialog: pick_folder, save_file

# ---------------------------------------------------------------------------
# Gui module split into src/Gui/ — see docs/gui.md for the file map.
# ---------------------------------------------------------------------------

include("Gui/BadAndStyling.jl")
include("Gui/State.jl")
include("Gui/PlotPanel.jl")
include("Gui/Persistence.jl")
include("Gui/TreePanel.jl")
include("Gui/InfoModal.jl")
include("Gui/Layout.jl")
