module Browser

include("Browser/MakieIntegration.jl")

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
