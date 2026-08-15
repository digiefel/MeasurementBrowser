# DataBrowser Architecture

## What this is

DataBrowser is a Julia working environment for data projects. It should be strictly
better than opening a REPL, finding files, loading them, extracting useful tables, computing values,
and writing figure code by hand. The app makes that workflow interactive: open a project, browse the
source structure, select collections or items, inspect data-derived values, and switch plots or views
without repeating the same parsing work.

Project code should describe the project in the same terms a script would use: which logical items a
source item contains, how to load data for those items, and how to present that data. The
package owns scanning, cache storage, background jobs, and UI state.

The intended experience is live and composable. A project can stay open while source files are added
or changed, and the browser should update without forcing the user back through startup or manual
reload steps. Views should be able to follow selections, collections, or matching rules, so a plot or
inspection tool can continue to show the relevant data as the source tree changes. Built-in
visualizers should handle common inspection tasks, while project code adds only the interpretation
and presentation details that are specific to the experiment.

## Package map

A solid arrow from A to B means A depends on B at compile time: A needs B to build.
A dotted arrow from A to B means A calls functions whose methods live in B through dispatch, even though A has no dependency on B.
🦆 marks DuckDB, 📈 GLMakie, 🖼 CImGui / GLFW.

```mermaid
flowchart TB
    db["<b>DataBrowser</b><br/>install target<br/>defines the public API"]
    prof["<b>DataBrowserProfiling</b><br/>internal tooling"]

    subgraph projects["Project packages"]
        recipes["<b>DataBrowserRecipes</b><br/>user-friendly register_*! API"]
        userpkg["⋯<br/>any user-authored<br/>project package"]
    end

    subgraph frontend["Frontend"]
        gui["<b>DataBrowserGUI</b> 🖼<br/>main GUI functionality"]
        plots["<b>DataBrowserPlots</b> 📈<br/>GLMakie-based plots"]
    end

    subgraph engine["Engine"]
        core["<b>DataBrowserCore</b><br/>workspace, job scheduling, coordination"]
        src["<b>DataBrowserSources</b><br/>DirectorySource, ..."]
        cache["<b>DataBrowserCache</b> 🦆<br/>Caching and persistence"]
        ann["<b>DataBrowserAnnotations</b><br/>tags, notes, layout"]
    end

    subgraph contracts["Contracts"]
        api["<b>DataBrowserAPI</b><br/>AbstractProject, stage and source interfaces, ItemIndex"]
    end

    db --> recipes
    db --> plots
    gui -->|"open_workspace<br/>modify_workspace!<br/>select_items!<br/>materialize_items<br/>workspace_status"| core
    gui --> cache
    gui --> ann
    plots --> gui
    plots -->|"read_item_data<br/>InspectorTable"| core
    core --> src
    core --> cache
    core --> ann
    core --> api
    src --> api
    cache --> api
    ann --> api
    recipes --> api
    userpkg --> api
    prof --> api

    core -.->|"read<br/>entries<br/>process<br/>analyze"| projects
    gui -.->|"draw!<br/>menu!<br/>init!"| plots
    cache -.->|"@timed_dbg"| prof
    core -.->|"@timed_dbg"| prof
    gui -.->|"@timed_dbg"| prof
    plots -.->|"@timed_dbg"| prof
    recipes -.->|"@timed_dbg"| prof

    classDef umbrellaC fill:#f7c9b8,stroke:#b5623f,color:#2b2b2b
    classDef plotsC fill:#f2c6d4,stroke:#b0466a,color:#2b2b2b
    classDef guiC fill:#d9c6ec,stroke:#7d54b0,color:#2b2b2b
    classDef coreC fill:#bcd4f0,stroke:#3f6fb0,color:#2b2b2b
    classDef recipesC fill:#f5ecc0,stroke:#b09a3f,color:#2b2b2b
    classDef srcC fill:#c9e4c5,stroke:#5a9a52,color:#2b2b2b
    classDef cacheC fill:#f5d6a8,stroke:#c08a3f,color:#2b2b2b
    classDef annC fill:#b8e0d2,stroke:#4f9a86,color:#2b2b2b
    classDef profC fill:#ded7ef,stroke:#8a7fb0,color:#2b2b2b
    classDef apiC fill:#cfd8e3,stroke:#5a6b80,color:#1e1e1e
    classDef ghostC fill:#f2f2f2,stroke:#9aa0a6,stroke-dasharray:5 4,color:#5f6368

    class userpkg ghostC
    class db umbrellaC
    class plots plotsC
    class gui guiC
    class core coreC
    class recipes recipesC
    class src srcC
    class cache cacheC
    class ann annC
    class prof profC
    class api apiC
```

## Core Flow

```
source item → interpret → logical data → process → analyze → collection process/analyze → views
                 │             │             │
                 └─ index      └─ DuckDB     └─ DuckDB + item metadata
```

A project/source implementation defines:

- interpreting each source item into logical data items
- processing one interpreted item
- computing per-item and per-collection metadata (item/collection `analyze`, collection `process`)
- defining project-specific visualizers when generic ones are not enough

The workspace owns:

- the open source(s) and their identity
- the progressively populated item index
- selection identities
- cache identity, freshness, storage, and repair
- scanning, cache work, progress, errors, and cancellation
- work dependency graph state and source fallback

The browser owns windows, controls, filters, and temporary rendering state. Annotations store
user-authored tags, notes, and other user-authored metadata. Other package modules own generic
visualizers, workflow persistence, and figure composition. User code should not know whether data
came from memory, cache, or the source. Package code does not know the meaning of a source item
beyond the contract methods it calls.

## Subpackages

Each package gets two views. A class diagram shows its data model: the structs, their fields, and
how they subtype the API abstractions and compose. The call flow diagram shows which functions
call which as arrows, with A -data-> B meaning "A calls B, which returns `data`".

A solid arrow is a call, labelled with the value the callee returns.
External calls from other packages are shown as thicker arrows.
"hot path" calls which are very frequent and performance-critical are shown in red.
A dashed arrow is data written to or read from storage.

A thick border marks a symbol the package exports.
A `?` after a field type means the field may be `nothing`.

### DataBrowserSources

#### Data model

```mermaid
%%{init: {'theme':'base','themeVariables':{'textColor':'#111','lineColor':'#555'}}}%%
classDiagram
    direction LR
    class AbstractDataSource:::api
    class AbstractDataSourceItem:::api
    class AbstractCollection:::api

    class DirectorySource:::srcExp {
        root_path : String
        recursive : Bool
        metadata_file : String?
        collection_metadata_entries : Dict
        has_metadata : Bool
        metadata_lock : ReentrantLock
        watcher_task : Task?
        watcher_cancel : CancellationTokenSource?
    }
    class SourceFile:::srcExp {
        filepath : String
        filename : String
        relative_path : String
        timestamp : DateTime?
        fingerprint : FileFingerprint
    }
    class FileFingerprint:::srcExp {
        path : String
        size_bytes : Int64
        mtime_ns : Int64
    }
    class DirectoryCollection:::srcInt {
        name : String
        metadata : Dict
    }

    AbstractDataSource <|-- DirectorySource
    AbstractDataSourceItem <|-- SourceFile
    AbstractCollection <|-- DirectoryCollection
    SourceFile *-- FileFingerprint : fingerprint

    classDef api fill:#cfd8e3,stroke:#5a6b80,stroke-width:1px,color:#111;
    classDef srcExp fill:#c9e4c5,stroke:#5a9a52,stroke-width:3px,color:#111;
    classDef srcInt fill:#c9e4c5,stroke:#5a9a52,stroke-width:1px,color:#111;
```

#### Call flow

```mermaid
%%{init: {'theme':'base','themeVariables':{'textColor':'#111','lineColor':'#555'}}}%%
flowchart TD
    subgraph sgScan["scan"]
        si("source_items")
        csf("collect_source_files")
        asf("append_source_files!")
        isf("index_source_file")
        ffp("file_fingerprint")
        isfn("is_source_filename")
        pts("parse_timestamp")
    end

    subgraph sgMeta["metadata.txt"]
        lcm("load_collection_metadata!")
        lcme("load_collection_metadata_entries")
        cmfp("collection_metadata_file_path")
        pmv("parse_metadata_value")
        ocm("own_collection_metadata")
        np("_named_path")
        dcp("default_collection_path")
        al("_annotated_level")
        acp("annotate_collection_path")
        MD[("collection_metadata_entries")]
    end

    subgraph sgLife["watch / lifecycle"]
        osrc("open_source")
        wsrc("watch_source")
        csrc("close_source!")
        cpyds("copy")
    end

    subgraph sgContract["contract methods"]
        idsf("id(::SourceFile) → String")
        lblsf("label(::SourceFile) → String")
        fpsf("fingerprint(::SourceFile) → FileFingerprint")
        sipf("source_item_path(::SourceFile) → String")
        sitf("source_item_timestamp(::SourceFile) → DateTime?")
        metasf("metadata(::SourceFile) → Dict")
        sid("source_id(::DirectorySource) → String")
        slbl("source_label(::DirectorySource) → String")
        snoun("source_item_noun(::DirectorySource) → String")
        iddc("id(::DirectoryCollection) → String")
        lbldc("label(::DirectoryCollection) → String")
        metadc("metadata(::DirectoryCollection) → Dict")
        recon("reconstruct(::DirectoryCollection) → DirectoryCollection")
    end

    extCore["DataBrowserCore"]
    extAPI["DataBrowserAPI"]
    extCache["DataBrowserCache"]
    extRecipes["DataBrowserRecipes"]
    extGUI["DataBrowserGUI"]
    extPlots["DataBrowserPlots"]

    subgraph sgExt["external"]
        FSYS["filesystem"]
        BFW["BetterFileWatching"]
        CT["CancellationTokens"]
    end

    si -->|"SourceFile[]"| csf
    csf --> asf
    asf -->|"SourceFile"| isf
    asf -->|"Bool"| isfn
    isf -->|"FileFingerprint"| ffp
    isf -->|"DateTime?"| pts
    csf -->|"String?"| cmfp
    csf -->|"walkdir"| FSYS
    ffp -->|"stat"| FSYS

    lcm -->|"String?"| cmfp
    lcm -->|"Dict"| lcme
    lcm -.->|"writes"| MD
    lcme -->|"Any"| pmv
    lcme -->|"readlines"| FSYS
    dcp -->|"AbstractCollection[]"| np
    np -->|"Dict"| ocm
    np -.->|"reads"| MD
    al -->|"Dict"| ocm
    al -.->|"reads"| MD
    al -->|"AbstractCollection"| recon
    al -->|"Dict"| metadc
    al -->|"String"| iddc
    acp -->|"AbstractCollection"| al
    acp -->|"String"| iddc

    osrc -->|"Bool"| lcm
    osrc -->|"isdir"| FSYS
    wsrc -->|"SourceFile[]"| csf
    wsrc -->|"Bool"| lcm
    wsrc -->|"String"| idsf
    wsrc -->|"FileFingerprint"| fpsf
    wsrc -->|"watch_folder"| BFW
    wsrc -->|"token"| CT
    csrc -->|"cancel"| CT

    extCore ==> osrc
    extCore ==> si
    extCore ==> wsrc
    extCore ==> csrc
    extCore ==> isf
    extCore ==> dcp
    extCore ==> acp
    extCore ==> sid
    extCore ==> snoun
    extCore ==> idsf
    extCore ==> lblsf
    extCore ==> metasf
    extCore ==> sipf
    extCore ==> sitf
    extCore ==> fpsf
    extCore ==> cpyds
    extCache ==> sid
    extCache ==> fpsf
    extCache ==> sipf
    extCache ==> sitf
    extAPI ==> sid
    extAPI ==> slbl
    extAPI ==> recon
    extRecipes ==> sipf
    extGUI ==> slbl
    extPlots ==> slbl

    wsrc -.->|"SourceChanges<br/>SourceError"| extCore

    idsf ~~~ lblsf ~~~ fpsf ~~~ sipf ~~~ sitf ~~~ metasf ~~~ sid ~~~ slbl ~~~ snoun ~~~ iddc ~~~ lbldc ~~~ metadc ~~~ recon ~~~ eqfp

    classDef srcExp fill:#c9e4c5,stroke:#5a9a52,stroke-width:3px,color:#111;
    classDef srcInt fill:#c9e4c5,stroke:#5a9a52,stroke-width:1px,color:#111;
    classDef ext fill:#eceff1,stroke:#90a4ae,stroke-width:1px,color:#111;
    classDef coreC fill:#bcd4f0,stroke:#3f6fb0,color:#111;
    classDef apiC fill:#cfd8e3,stroke:#5a6b80,color:#111;
    classDef cacheC fill:#f5d6a8,stroke:#c08a3f,color:#111;
    classDef recipesC fill:#f5ecc0,stroke:#b09a3f,color:#111;
    classDef guiC fill:#d9c6ec,stroke:#7d54b0,color:#111;
    classDef plotsC fill:#f2c6d4,stroke:#b0466a,color:#111;

    class isf,ffp srcExp;
    class si,csf,asf,isfn,pts,lcm,lcme,cmfp,pmv,ocm,np,dcp,al,acp,MD,osrc,wsrc,csrc,cpyds,idsf,lblsf,fpsf,sipf,sitf,metasf,sid,slbl,snoun,iddc,lbldc,metadc,recon,eqfp srcInt;
    class extCore coreC;
    class extAPI apiC;
    class extCache cacheC;
    class extRecipes recipesC;
    class extGUI guiC;
    class extPlots plotsC;
    class FSYS,BFW,CT ext;

    linkStyle 2,4,8 stroke:#d1495b,stroke-width:3px;
```
