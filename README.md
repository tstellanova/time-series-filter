# time-series-filter

Some convenient traits and types for
processing sequential data.

- `FloatSeriesEwmaFilter` can be used to track the 
exponential weighted moving average (EWMA) of a varying
signal.  This is essentially an infinite impulse response (IIR)
filter and a low pass filter (LPF).
- `IntSeriesEwmaFilter` creates a filter for integer types, 
which avoids floating point math.


## Examples

See the tests in [lib.rs](./src/lib.rs) for examples of usage. 


## License

BSD-3-Clause, see `LICENSE` file. 