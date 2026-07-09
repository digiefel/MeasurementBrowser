import CImGui as ig

"""Render a small question mark with explanatory hover text."""
function _helpmarker(description::String)::Nothing
    ig.TextDisabled("(?)")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(description)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
    return nothing
end
