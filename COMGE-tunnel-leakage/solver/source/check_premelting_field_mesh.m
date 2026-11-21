function evidence = check_premelting_field_mesh(model)
%CHECK_PREMELTING_FIELD_MESH Necessary trace-resolution check before UQ.
% Full rank is NOT an inf-sup proof or a substitute for grid convergence.
radialCells=numel(unique(model.mesh.nodeRingIndex))-1;
B=full(model.gapOperator); traceRank=rank(B,1e-12);
evidence=struct('mechanicalRadialCells',radialCells, ...
    'hydraulicCellsPerJoint',model.nq,'gapTraceRank',traceRank,'traceRows',size(B,1), ...
    'radialResolutionAdequate',radialCells>=model.nq,'traceFullRowRank',traceRank==size(B,1));
assert(evidence.radialResolutionAdequate && evidence.traceFullRowRank, ...
    'PhaseSwitch:TraceResolution','Refine the segment radial mesh together with the fluid cells; pressure cells exceed the resolved opening trace.');
end
