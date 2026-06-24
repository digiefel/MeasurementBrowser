module Browser

include("Browser/MakieIntegration.jl")

include("Browser/State.jl")
include("Browser/Performance.jl")
include("Browser/Tags.jl")
include("Browser/Persistence.jl")
include("Browser/Operations.jl")
include("Gui/Utilities.jl")
include("Gui/TagStyling.jl")
include("Gui/PlotPanel.jl")
include("Gui/DataGrid.jl")
include("Gui/TableInspector.jl")
include("Gui/TreePanel.jl")
include("Gui/InfoModal.jl")
include("Gui/Layout.jl")
include("Gui/PerformanceWindow.jl")

export open_browser

end
