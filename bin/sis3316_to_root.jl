#!/usr/bin/env julia

# This file is a part of SIS3316.jl, licensed under the MIT License (MIT).


using SIS3316
using ROOTFramework


mktemp_custom(parent=tempdir(), format="tmpXXXXXX") = begin
    b = joinpath(parent, format)
    p = ccall(:mkstemp, Int32, (Ptr{UInt8},), b) # modifies b
    systemerror(:mktemp, p == -1)
    return (b, fdio(p, true))
end


sis3316_to_root(input_io::IO, output_tfile::TFile; evt_merge_window::AbstractFloat = 100e-9) = begin
    Base.time(x::Pair{Int64, SIS3316.RawChEvent}) = time(x.second)

    const ttoutput = TTreeOutput("data", "Data")

    const info_idx = ttoutput[:info_idx] = Ref{Int32}(0)
    const info_time = ttoutput[:info_time] = Ref{Float64}(0)

    const raw_pp_ch = ttoutput[:raw_pp_ch] = Vector{Int32}()
    const raw_pp_mca = ttoutput[:raw_pp_mca] = Vector{Int32}()
    const raw_pp_trig_max = ttoutput[:raw_pp_trig_max] = Vector{Int32}()
    const raw_pp_peak_pos = ttoutput[:raw_pp_peak_pos] = Vector{Int32}()
    const raw_pp_peak_height = ttoutput[:raw_pp_peak_height] = Vector{Int32}()
    const raw_pp_acc = map(i -> ttoutput[Symbol("raw_pp_acc_$i")] = Vector{Int32}(), 1:8)

    const raw_trig_ch = ttoutput[:raw_trig_ch] = Vector{Int32}()
    const raw_trig_trel = ttoutput[:raw_trig_trel] = Vector{Float64}()
    const raw_trig_pileup = ttoutput[:raw_trig_pileup] = Vector{Int32}()
    const raw_trig_overflow = ttoutput[:raw_trig_overflow] = Vector{Int32}()

    const ch_sized_vecs = Vector[raw_pp_ch, raw_pp_mca, raw_pp_trig_max,
        raw_pp_peak_pos, raw_pp_peak_height,
        raw_trig_ch, raw_trig_trel, raw_trig_pileup, raw_trig_overflow
    ]
    append!(ch_sized_vecs, raw_pp_acc)

    for v in ch_sized_vecs sizehint!(v, 16) end


    const reader = eachchunk(input_io, SIS3316.UnsortedEvents)
    open(ttoutput, output_tfile)

    local evtno = 0

    for unsorted in reader
        const sorted = sortevents(unsorted, merge_window = evt_merge_window)
        const evtv = Vector{Pair{Int64, SIS3316.RawChEvent}}()
        const timestamps = Vector{Float64}()

        const energynull = SIS3316.EnergyValues(0, 0)
        const mawnull = SIS3316.MAWValues(0, 0, 0)
        const psanull = SIS3316.PSAValue(0, 0)
        const flagsnull = SIS3316.EvtFlags(false,false,false,false)

        for evt in sorted
            evtno += 1

            resize!(evtv, length(evt))
            copy!(evtv, evt)
            sort!(evtv, by = first)
            resize!(timestamps, length(evtv))
            map!(time, timestamps, evtv)
            const starttime = isempty(timestamps) ? zero(Float64) : minimum(timestamps)

            info_idx.x = evtno
            info_time.x = starttime

            for v in ch_sized_vecs empty!(v) end

            for (ch, chevt) in evtv
                push!(raw_pp_ch, ch)
                push!(raw_pp_mca, get(chevt.energy, energynull).maximum)
                push!(raw_pp_trig_max, get(chevt.trig_maw, mawnull).maximum)
                push!(raw_pp_peak_pos, get(chevt.peak_height, psanull).index)
                push!(raw_pp_peak_height, get(chevt.peak_height, psanull).value)

                for i in eachindex(raw_pp_acc)
                    push!(raw_pp_acc[i], get(chevt.accsums, i, 0))
                end

                push!(raw_trig_ch, ch)
                push!(raw_trig_trel, time(chevt) - starttime)
                push!(raw_trig_pileup, chevt.pileup_flag + 2 * get(chevt.flags, flagsnull).pileup +  4 * get(chevt.flags, flagsnull).repileup)
                push!(raw_trig_overflow, 1 * get(chevt.flags, flagsnull).overflow +  2 * get(chevt.flags, flagsnull).underflow)
            end

            push!(ttoutput)
        end
    end
end


sis3316_to_root(input_fname::AbstractString; evt_merge_window::AbstractFloat = 100e-9) = begin
    const fnexpr = r"(.*)\.dat(\.[^.]+)?"
    const fnbase = match(fnexpr, basename(input_fname))[1]
    const output_fname = "$(fnbase).root"
    output_tmpname, tmpio = mktemp_custom(pwd(), "$(output_fname).tmp-XXXXXX")
    close(tmpio)

    const input_io = open_decompressed(input_fname)
    const output_tfile = TFile(output_tmpname, "recreate")

    sis3316_to_root(input_io, output_tfile, evt_merge_window = evt_merge_window)

    close(input_io)
    close(output_tfile)

    mv(output_tmpname, output_fname, remove_destination = true)
end



main() = begin
    local evt_merge_window = 100e-9

    const inputs = ARGS

    for input_fname in inputs
        try
            info("Converting \"$(input_fname)\"")
            @time sis3316_to_root(input_fname, evt_merge_window = evt_merge_window)
        catch err
            print_with_color(:red, STDERR, "ERROR: $err\n")
        end
    end
end

main()
