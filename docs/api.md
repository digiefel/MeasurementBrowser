# Public project API

The user-facing project documentation is organized as a book under [`docs/public/`](public/index.md).

- [How DataBrowser works](public/pipeline.md)
- [Your first project](public/getting-started.md)
- [Registration API](public/registration.md)
- [Metadata and collections](public/metadata-and-collections.md)
- [Type API](public/type-api.md)
- [API reference](public/reference.md)

The book is built with Documenter by running:

```bash
julia --project=docs docs/make.jl
```
