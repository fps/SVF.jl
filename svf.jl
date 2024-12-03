module StateVariableFilter

struct CommonProcessingResult{T<:AbstractFloat}
  w::T
  a::T
  b::T
  c1::T
  c2::T
end

function process_common(freq, q)
  w = 2.0 * tan(pi * freq)
  a = w / q
  b = w*w
  c1 = (a + b) / (1 + 0.5 * a + 0.25 * b)
  c2 = b / (a + b)

  CommonProcessingResult(w, a, b, c1, c2)
end

mutable struct State{T<:AbstractFloat}
  z1::T
  z2::T
end

State() = State(0.0,0.0)

struct D1D0{T<:AbstractFloat}
  d1::T
  d0::T
end

function calculate_d0_high(common)
  1.0 - 0.5 * common.c1 + common.c1 * common.c2 * 0.25
end

function calculate_d1d0_band(common)
  d1_band = 1.0 - common.c2
  d0_band = d1_band * common.c1 * 0.5
  D1D0(d1_band, d0_band)
end

function calculate_d0_low(common)
  common.c1 * common.c2 * 0.25
end

function process(s::State, in, highgain, bandgain, lowgain, freq, q)
    common = process_common(freq, q)

    d0_high = calculate_d0_high(common)
    d1d0_band = calculate_d1d0_band(common)
    d0_low = calculate_d0_low(common)

    out = 0

    x = in - s.z1 - s.z2 + 1e-20
    out += highgain * d0_high * x
    out += bandgain * (d1d0_band.d0 * x + d1d0_band.d1 * s.z1)
    s.z2 += common.c2 * s.z1
    out += lowgain * (d0_low * x + s.z2)
    s.z1 += common.c1 * x

    out
end


function process!(s::State, input_buffer, output_buffer, sample_count, highgain, bandgain, lowgain, freq, q)
    common = process_common(freq, q)

    d0_high = calculate_d0_high(common)
    d1d0_band = calculate_d1d0_band(common)
    d0_low = calculate_d0_low(common)
    
    for sample_index = 1:sample_count
        in = input_buffer[sample_index]
        out = 0

        x = in - s.z1 - s.z2 + 1e-20
        out += highgain * d0_high * x
        out += bandgain * (d1d0_band.d0 * x + d1d0_band.d1 * s.z1)
        s.z2 += common.c2 * s.z1
        out += lowgain * (d0_low * x + s.z2)
        s.z1 += common.c1 * x

        output_buffer[sample_index] = out;
    end
    
    output_buffer
end
