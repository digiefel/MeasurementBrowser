"""CImGui browser shell and panels."""
module DataBrowserGUI

include("Browser.jl")

using .Browser: open_browser, close_browser!, BrowserSession
export open_browser, close_browser!, BrowserSession, Browser

end
