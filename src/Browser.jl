module Browser

include("Browser/MakieIntegration.jl")

include("Browser/State.jl")
include("Browser/Performance.jl")
include("Gui/BadAndStyling.jl")
include("Browser/FigureScripts.jl")
include("Browser/Persistence.jl")
include("Browser/Operations.jl")
include("Gui/Utilities.jl")
include("Gui/PlotPanel.jl")
include("Gui/TreePanel.jl")
include("Gui/InfoModal.jl")
include("Gui/Layout.jl")

export start_browser

end
