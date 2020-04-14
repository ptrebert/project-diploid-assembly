
Queue: default
    queue_type = Route
    total_jobs = 0
    state_count = Transit:0 Queued:0 Held:0 Waiting:0 Running:0 Exiting:0 Begun
        :0 
    resources_max.walltime = 167:59:59
    route_destinations = CUDA,short,workq,long
    route_retry_time = 30
    enabled = True
    started = True


Queue: short
    queue_type = Execution
    Priority = 120
    total_jobs = -14
    state_count = Transit:0 Queued:0 Held:1 Waiting:0 Running:4 Exiting:0 Begun
        :0 
    max_queued = [u:PBS_GENERIC=1024]
    from_route_only = True
    resources_max.ngpus = 0
    resources_max.walltime = 02:00:00
    resources_default.preempt_targets = NONE
    resources_assigned.mem = 11534336kb
    resources_assigned.mpiprocs = 0
    resources_assigned.ncpus = 4
    resources_assigned.nodect = 4
    max_run = [u:PBS_GENERIC=512]
    enabled = True
    started = True


Queue: workq
    queue_type = Execution
    Priority = 100
    total_jobs = 182
    state_count = Transit:0 Queued:0 Held:11 Waiting:0 Running:105 Exiting:0 Be
        gun:1 
    max_queued = [u:PBS_GENERIC=1024]
    from_route_only = True
    resources_max.ngpus = 0
    resources_max.walltime = 72:00:00
    resources_min.walltime = 02:00:01
    resources_default.preempt_targets = NONE
    resources_assigned.mem = 4672716800kb
    resources_assigned.mpiprocs = 2396
    resources_assigned.ncpus = 2626
    resources_assigned.nodect = 141
    enabled = True
    started = True


Queue: long
    queue_type = Execution
    Priority = 80
    total_jobs = 82
    state_count = Transit:0 Queued:26 Held:2 Waiting:0 Running:95 Exiting:0 Beg
        un:4 
    max_queued = [u:PBS_GENERIC=1024]
    from_route_only = True
    resources_max.ngpus = 0
    resources_max.walltime = 167:00:00
    resources_min.walltime = 72:00:01
    resources_default.preempt_targets = NONE
    resources_assigned.mem = 1525022720kb
    resources_assigned.mpiprocs = 24
    resources_assigned.ncpus = 272
    resources_assigned.nodect = 100
    max_run_res.ncpus = [o:PBS_ALL=1024]
    enabled = True
    started = True