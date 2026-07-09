"""
Concrete data sources: directory scanning and delimited-file preview.
"""
module DataBrowserSources

using DataBrowserAPI

include("directory_source.jl")
include("tabular_file_source.jl")

export DirectorySource, SourceFile, FileFingerprint
export DEFAULT_DIRECTORY_METADATA_FILE
export TabularFileSource, inspect_table
export index_source_file, file_fingerprint

end
