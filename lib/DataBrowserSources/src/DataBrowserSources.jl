"""
Concrete data sources: directory scanning.
"""
module DataBrowserSources

using DataBrowserAPI

include("directory_source.jl")

export DirectorySource, SourceFile, FileFingerprint
export DEFAULT_DIRECTORY_METADATA_FILE
export index_source_file, file_fingerprint

end
