using Documenter
using DocumenterMermaid
using Literate

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const GENERATED_SOURCE = joinpath(@__DIR__, "generated")
const GENERATED_CODE = joinpath(GENERATED_SOURCE, "code")
const EXCLUDED_SOURCE_PREFIXES = (".agents/", "bench/")
const EXCLUDED_SOURCE_FILES = Set(["docs/make.jl"])

function include_source(source_path::AbstractString)
    source_path in EXCLUDED_SOURCE_FILES && return false
    return !any(prefix -> startswith(source_path, prefix), EXCLUDED_SOURCE_PREFIXES)
end

function tracked_julia_files()
    command = Cmd(
        Cmd([
            "git", "ls-files", "--cached", "--others", "--exclude-standard", "--",
            "*.jl",
        ]);
        dir=REPOSITORY_ROOT,
    )
    return sort(filter(include_source, filter(!isempty, readlines(command))))
end

function page_path(source_path::AbstractString)
    stem, _ = splitext(source_path)
    return joinpath("code", stem * ".md")
end

function generate_source_page(source_path::AbstractString)
    input_path = joinpath(REPOSITORY_ROOT, source_path)
    output_path = joinpath(GENERATED_SOURCE, page_path(source_path))
    output_directory = dirname(output_path)
    mkpath(output_directory)

    title = "# # `$(source_path)`\n#\n"
    Literate.markdown(
        input_path,
        output_directory;
        name=splitext(basename(source_path))[1],
        flavor=Literate.CommonMarkFlavor(),
        execute=false,
        credit=false,
        preprocess=source -> title * source,
    )
    return page_path(source_path)
end

function copy_architecture_page()
    cp(
        joinpath(@__DIR__, "ARCHITECTURE.md"),
        joinpath(GENERATED_SOURCE, "ARCHITECTURE.md");
        force=true,
    )
end

function navigation_entry(source_page, prefix)
    source_path, markdown_path = source_page
    label = String(last(split(source_path, prefix; limit=2)))
    return label => markdown_path
end

function grouped_navigation(source_pages)
    root_package = [
        navigation_entry(page, "src/")
        for page in source_pages if startswith(first(page), "src/")
    ]

    library_names = sort(unique(String(split(first(page), '/')[2]) for page in source_pages if startswith(first(page), "lib/")))
    libraries = Any[
        library_name => [
            navigation_entry(page, "lib/$(library_name)/")
            for page in source_pages if startswith(first(page), "lib/$(library_name)/")
        ]
        for library_name in library_names
    ]

    example_names = sort(unique(String(split(first(page), '/')[2]) for page in source_pages if startswith(first(page), "examples/")))
    examples = Any[
        example_name => [
            navigation_entry(page, "examples/$(example_name)/")
            for page in source_pages if startswith(first(page), "examples/$(example_name)/")
        ]
        for example_name in example_names
    ]

    tests = [
        navigation_entry(page, "test/")
        for page in source_pages if startswith(first(page), "test/")
    ]

    return Any[
        "DataBrowser" => root_package,
        "Libraries" => libraries,
        "Examples" => examples,
        "Tests" => tests,
    ]
end

function write_navigation(io, entries, heading_level)
    for (label, destination) in entries
        if destination isa AbstractString
            markdown_path = replace(destination, '\\' => '/')
            println(io, "- [`$(label)`]($(markdown_path))")
        else
            println(io)
            println(io, repeat("#", heading_level), " ", label)
            write_navigation(io, destination, heading_level + 1)
        end
    end
end

function write_index(navigation)
    open(joinpath(GENERATED_SOURCE, "index.md"), "w") do io
        println(io, "# DataBrowser source")
        println(io)
        println(io, "[Architecture](ARCHITECTURE.md)")
        write_navigation(io, navigation, 2)
    end
end

rm(GENERATED_SOURCE; recursive=true, force=true)
mkpath(GENERATED_CODE)

source_paths = tracked_julia_files()
source_pages = [source_path => generate_source_page(source_path) for source_path in source_paths]
navigation = grouped_navigation(source_pages)
copy_architecture_page()
write_index(navigation)

makedocs(
    root=@__DIR__,
    source="generated",
    build="build",
    clean=true,
    sitename="DataBrowser source",
    format=Documenter.HTML(prettyurls=false),
    doctest=false,
    checkdocs=:none,
    pages=[
        "Source" => "index.md",
        "Architecture" => "ARCHITECTURE.md",
        navigation...,
    ],
)
