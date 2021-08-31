

rule check_overlong_edges:
    input:
        gfa = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.tigs]
    output:
        discard = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.discard.links'
    conda: '../../../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('gfa_check_ovl.py'),
        stats_out = lambda wildcards, output: output.discard.replace('.discard.links', '.stats')
    resources:
        mem_total_mb = lambda wildcards, attempt: 512 * attempt,
        runtime_hrs = 0,
        runtime_min = lambda wildcards, attempt: 10 * attempt,
    shell:
        '{params.script_exec} -g {input.gfa} > {params.stats_out}'


rule clean_input_gfa:
    input:
        gfa = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.tigs],
        discard = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.discard.links'
    output:
        gfa = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.gfa'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import os
        import io
        import shutil

        skip_lines = []
        with open(input.discard, 'r') as table:
            for line in table:
                if line.startswith('#') or not line.strip():
                    continue
                ln = line.split('\t')[0]
                skip_lines.append(int(ln))
        
        if not skip_lines:
            shutil.copy(input.gfa, output.gfa)
        else:
            gfa_buffer = io.StringIO()
            with open(input.gfa, 'r') as gfa:
                for ln, line in enumerate(gfa, start=1):
                    if ln in skip_lines:
                        continue
                    gfa_buffer.write(line)
            with open(output.gfa, 'w') as gfa:
                _ = gfa.write(gfa_buffer.getvalue())
    # END OF RUN BLOCK