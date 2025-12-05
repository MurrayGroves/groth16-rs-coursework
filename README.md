# Groth16

This crate provides a library for generating and verifying Groth16 proofs over a given Rank 1 Constraint System.

The implementation is type generic over any curve that implements `ark_ec::Pairing`, but is tested in this library with
MNT6-753 since it's the highest security MNT curve.

The polynomial implementation (including Lagrangian interpolation, SRS evaluation, scalar and polynomial addition,
subtraction,
multiplication, and division) was developed for this library due to my dissatisfaction with the Lagrangian interpolation
available in existing polynomial libraries, and is generic over any type that implements `ark_ff::Field`, with the SRS
evaluation being generic over any type implementing `ark_ec::CurveGroup`.

The public interface is well-documented and may be viewed at `docs/groth16/index.html`, or by running
`cargo doc --open`.

Additionally, errors are handled using the `rootcause` crate which allows for developer friendly backtraces including
additional context and attached debug information (think showing operands to a failed call). All Results are annotated
with at the minimum a `.context()` call describing the operation that failed.

The core types are `circuits::R1CS`, `circuits::QAP`, `groth16::TrustedSetupOutput`, `groth16::Proof`.

The general flow is:

- Define R1CS
- Use `QAP::from` to convert R1CS to QAP
- Generate a Trusted Setup using `TrustedSetupOutput::new`
- Generate a proof using `trusted_setup.prove(witness)`
- Verify proof with `proof.verify()`

Tests can be run on this library by running `cargo test`.

The library is extensively tested with 100% test coverage, which can be verified with
`cargo install tarpaulin; cargo tarpaulin --engine llvm` or by looking at the coverage report in `coverage`. Note that
although the coverage report reports 414/418 lines covered, the missing 4 lines are all match arms that are incorrectly
marked as being code, which can be verified by going into `polynomial.rs` in the report and finding the red lines.
