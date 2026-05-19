# Measurement Metadata Contract

## Goal

Make measurement metadata simple enough that users can trust it and developers can extend it without
guessing.

Every value shown by the browser must answer three questions clearly:

```text
What object does this describe?
Where did it come from?
Was it parsed directly or computed later?
```

This document defines the user-facing contract. It does not prescribe the internal implementation.
Here, "user" includes both people using the browser and people adding a new project or measurement
procedure.

## Buckets

Measurement metadata has three public buckets.

```text
Device parameters
  Facts about the physical device.
  Example: area, film thickness, coordinates.

Measurement procedure parameters
  Facts about how this measurement entry was configured or identified.
  Parsed cheaply from the filename, CSV header, or the local file structure.

Analyzed stats
  Facts computed from loaded data, cached analysis, or measurement history.
```

The same concept must not appear in multiple buckets with competing meanings.

## User Contract

- Device parameters describe the device, not the measurement.
- Measurement procedure parameters describe the measurement entry before analysis.
- Analyzed stats describe what the data or measurement history shows after analysis.
- Project/procedure authors must have one obvious place to declare the fields their procedure
  exposes.
- Project/procedure authors must be able to tell whether a field should be parsed as a procedure
  parameter or computed as an analyzed stat without reading UI, cache, or plotting code.
- Missing numeric values are explicit and consistent, not silently absent because one file type took
  a different code path.
- Adding an ordinary field should require one obvious declaration and one obvious test.
- UI and cache should preserve the bucket a value belongs to.
- Tests are written before implementation and define the intended shape.

## First Contract: PUND Family

PUND, wakeup, and fatigue are the first family to migrate because they currently mix filenames,
headers, selected entries, and computed values.

For this family, the browser must expose:

```text
Measurement procedure parameters
  pulse history and settings before the readout waveform

Analyzed stats
  actual values computed from the selected waveform or measurement history
```

The first required procedure parameters for every PUND-family logical measurement are:

```text
wakeup_count
wakeup_f
wakeup_V
fatigue_count
fatigue_f
fatigue_V
```

The count fields describe pulses already applied to the device before the readout waveform,
including pulses sent earlier in the same file. Count fields default to `0` when that kind of pulse
history is absent. Frequency and voltage fields default to `NaN` when they are unavailable.

Wakeup and fatigue files have different structures and must have separate extractors. They share
the public field names, not a generic parsing path. Wakeup counts accumulate per device and wakeup
condition across files. Fatigue measurements inherit prior wakeup history for the same device.

The first required analyzed voltage stats for every PUND-family logical measurement are:

```text
V_base
V_min
V_max
V_amp
```

These are stats because they describe the loaded waveform, not merely what the measurement was
requested to do.

The existing fixture tests are the contract for the first migration.

## Non-Goals

- Do not keep compatibility aliases for known-wrong PUND metadata names.
- Do not move device metadata into measurement parsing.
- Do not treat waveform-derived values as parsed measurement parameters.
- Do not make UI text assertions the source of truth for metadata correctness.

## Rollout

1. Add contract tests for the target metadata shape using real fixture files. Done.
2. Add the metadata storage fields to the data model. Done.
3. Move cheap local extraction into measurement procedure parameters.
4. Move waveform/history-derived values into analyzed stats.
5. Migrate PUND/wakeup/fatigue first and delete old mixed keys.
6. Use the same pattern for other procedure families later.

## Open Question

How much of this should align now with the future project/procedure API rewrite?
