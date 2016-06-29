#!/usr/bin/env julia

# This file is a part of SIS3316.jl, licensed under the MIT License (MIT).


using ArgParse
using SIS3316
using ROOTFramework


mktemp_custom(parent=tempdir(), format="tmpXXXXXX") = begin
    b = joinpath(parent, format)
    p = ccall(:mkstemp, Int32, (Ptr{UInt8},), b) # modifies b
    systemerror(:mktemp, p == -1)
    return (b, fdio(p, true))
end


sis3316_to_root(input_io::IO, output_tfile::TFile; evt_merge_window::AbstractFloat = 100e-9, waveform_format = :none) = begin
    Base.time(x::Pair{Int64, SIS3316.RawChEvent}) = time(x.second)

    const ttoutput = TTreeOutput("data", "Data")

    const info_idx = ttoutput[:info_idx] = Ref{Int32}(0)
    const info_time = ttoutput[:info_time] = Ref{Float64}(0)

    const raw_pp_ch = ttoutput[:raw_pp_ch] = Int32[]
    const raw_pp_mca = ttoutput[:raw_pp_mca] = Int32[]
    const raw_pp_trig_max = ttoutput[:raw_pp_trig_max] = Int32[]
    const raw_pp_peak_pos = ttoutput[:raw_pp_peak_pos] = Int32[]
    const raw_pp_peak_height = ttoutput[:raw_pp_peak_height] = Int32[]
    const raw_pp_acc = map(i -> ttoutput[Symbol("raw_pp_acc_$i")] = Int32[], 1:8)

    const raw_trig_ch = ttoutput[:raw_trig_ch] = Int32[]
    const raw_trig_trel = ttoutput[:raw_trig_trel] = Float64[]
    const raw_trig_pileup = ttoutput[:raw_trig_pileup] = Int32[]
    const raw_trig_overflow = ttoutput[:raw_trig_overflow] = Int32[]

    const raw_wf_ch = Int32[]
    const raw_wf_smpl_n = Int32[]
    const raw_wf_smpl_v = Int32[]

    if waveform_format == :none
    elseif waveform_format == :integers
        ttoutput[:raw_wf_ch] = raw_wf_ch
        ttoutput[:raw_wf_smpl_n] = raw_wf_smpl_n
        ttoutput[:raw_wf_smpl_v] = raw_wf_smpl_v
    else
        error("Unknown waveform format $waveform_format.")
    end

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

            empty!(raw_wf_ch)
            empty!(raw_wf_smpl_n)
            empty!(raw_wf_smpl_v)

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

                if !isempty(chevt.samples)
                    if waveform_format == :integers
                        push!(raw_wf_ch, ch)
                        push!(raw_wf_smpl_n, length(chevt.samples))
                        append!(raw_wf_smpl_v, chevt.samples)
                    end
                end
            end

            push!(ttoutput)
        end
    end
end


sis3316_to_root(input_fname::AbstractString; evt_merge_window::AbstractFloat = 100e-9, waveform_format = :none) = begin
    const fnexpr = r"(.*)\.dat(\.[^.]+)?"
    const fnbase = match(fnexpr, basename(input_fname))[1]
    const output_fname = "$(fnbase).root"

    if !ispath(output_fname)
        output_tmpname, tmpio = mktemp_custom(pwd(), "$(output_fname).tmp-XXXXXX")
        close(tmpio)

        const input_io = open_decompressed(input_fname)
        const output_tfile = TFile(output_tmpname, "recreate")

        sis3316_to_root(input_io, output_tfile, evt_merge_window = evt_merge_window, waveform_format = waveform_format)

        close(input_io)
        close(output_tfile)

        mv(output_tmpname, output_fname, remove_destination = false)
    else
        info("Output file \"$(output_fname)\" already exists, skipping \"$(input_fname)\".")
    end
end



main() = begin
    argsettings = ArgParseSettings()

    @add_arg_table argsettings begin
        "--merge-within", "-m"
            help = "time window for merging events (no merging for negative values)"
            arg_type = Float64
            default = 100e-9
        "--waveform-format", "-w"
            help = "waveform format (\"none\" or \"integers\")"
            arg_type = AbstractString
            default = "integers"
        "inputs"
            help = "input files"
            nargs = '*'
    end

    const parsed_args = parse_args(argsettings)

    for input_fname in parsed_args["inputs"]
        try
            info("Converting \"$(input_fname)\"")
            @time sis3316_to_root(
                input_fname,
                evt_merge_window = parsed_args["merge-within"],
                waveform_format = symbol(parsed_args["waveform-format"])
            )
        catch err
            print_with_color(:red, STDERR, "ERROR: $err\n")
        end
    end
end

main()
