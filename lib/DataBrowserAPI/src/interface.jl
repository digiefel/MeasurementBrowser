"""Classify the item kind represented by a project source filename."""
function detect_kind end

"""Return the human-readable label for a project item kind."""
function kind_label end

"""Return the human-readable label for one logical item."""
function display_label end

"""Return the project-specific display label for one collection path."""
function collection_path_label end
