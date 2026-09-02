# Guidelines

IMPORTANT ABOVE ALL ELSE: communicate following ISO 24495-1, i.e. plain language practices.
Plain language is communication that puts readers first. It considers:
— what readers want and need to know;
— readers’ level of interest, expertise and literacy skills;
— **the context in which readers will use the document**.
Plain language ensures readers can find what they need, understand it and use it. Thus, plain language focuses on how successfully readers can use the document rather than on mechanical measures such as readability formulas.
Extensive studies have shown that writing in plain language saves time or money (or both) for readers and organizations. Plain language is more effective and produces better outcomes. In addition, readers prefer plain language. For organizations, plain language is an important way to build trust with the readers. Finally, the process of translating is more efficient for plain language documents than for documents that are difficult to understand.
Plain language is not to be confused with easy language. Plain language can be used for a general audience, while easy language is used for people who have difficulties with reading comprehension. 

Do not use meta-language or "punchy" figures of speech. Do not say that you will consider the context in which readers will use the document, or "announce" your intent. Think about the intent, think about your output, but produce actual prose/code without self-commentary or meta-references.

Keep consistent vocabulary. Avoid synonims and colorful prose. Always be specific in your wording.
Avoid overusing internal jargon. Reduce cognitive load by explaining jargon when useful.

IMPORTANT: before planning multi-package changes, and whenever broad context is needed, read the north-star document: [docs/vision.md](docs/vision.md).
For the full architectural model, when needed, see [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).
When making a change or looking to add a feature, read the roadmap: [docs/roadmap.md](docs/roadmap.md).
Benchmark details: [bench/README.md](bench/README.md).

Use docstrings when useful, and ALWAYS have docstrings on public APIs. 

This is a pre-alpha with zero users: no compatibility or migration code is ever needed.
Refactors are encouraged whenever tension arises, as few models are fixed in stone.
When code and docs disagree, fix the doc in the same commit.

## Commands

```bash
# Run tests and write bench/status.txt (skip for doc-only / trivial edits; time consuming)
julia --project=bench --threads=auto test/runtests.jl

# Generate public docs
julia --project=docs docs/make.jl
```

## Architecture
Project scripts describe how to recognize files, parse them into items, and draw plots. The package  
handles directory scanning, background processing, DuckDB caching, the item tree, selection, and the  
browser UI.  
Project code should not touch cache files, background jobs, or UI state. It is user code: a black box,  
never to be thought about or optimized.

When a workspace opens, the package works through a dependency graph with multiple stages: 
interpret source items, process each item, analyze each item, then process and analyze at the collection level.
If the user selects items that are still processing, that work gets higher priority. 
Full detail: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

Planned work lives in `docs/roadmap.md`.

## Testing

Before a commit, run the full suite once:
`julia --project=bench --threads=auto test/runtests.jl` (or `bench/run.sh` if you want peak RSS in
`status.txt`). That command uses the `bench/` environment for unit tests and then writes
`bench/status.txt`. Skip for doc-only, inspection-only, or harmless local edits. Fixtures in
`test/fixtures/`; the inline project lives in `test/test_project.jl`. Plot/GUI tests: metadata,
labels, figure creation — not pixels.

## Benchmarks

The performance run is the last step of the test command above. See [bench/README.md](bench/README.md).

## Work style

Work in small reviewable items. Always clear up confusion. Do not assume. Ask the user whenever there's decisions,
proposing various options to stimulate ideas. Asking is always better than assuming.
Conversations are always preferred to long outputs. A question mark is worth 1000 words.

After a turn with changes still pending, propose a series of one or more commits with a simple title.
Remember to advise the user when a commit or more are due, or overdue.
