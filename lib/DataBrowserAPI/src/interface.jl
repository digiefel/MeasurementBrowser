"""Return the stable name used to identify a project."""
function project_name end

"""Return a short human-readable description of a project."""
function project_description end

"""Classify the item kind represented by a project source filename."""
function detect_kind end

"""Return the human-readable label for a project item kind."""
function kind_label end

"""Return the human-readable label for one logical item."""
function display_label end

"""Return the project-specific display label for one collection path."""
function collection_path_label end

"""Clear any per-scan timing a project accumulates. Called once at the start of every scan."""
function reset_scan_profile! end

"""Record one callback phase for a source item in the current scan."""
function record_scan_phase! end

"""Finish one source-item timing after all expanded items have been processed."""
function finish_source_profile! end

"""
Per-item-kind timing for the most recent scan, newest scan replacing the last.

Returns one aggregate row per item kind. Surfaced in the performance window.
"""
function scan_profile_summary end

"""Return source-item timing rows for the latest scan, slowest first."""
function scan_source_profile end
