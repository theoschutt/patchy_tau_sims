for i in {200..980..20}
do
    let version=$i/20
    printf -v v_ext "%03d" $version
    echo "Starting mock batch: $version"
    python run_mocks_pipeline.py --nmocks 20 --njobs 4 --start_seed $i --version $v_ext --equalsignedweights > hi-tsz+lens_final_v2_20mocks_${v_ext}.log
    echo "Mock batch $version done."
    echo "Deleting batch flatmaps..."
    trash /home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_final_v2_20mocks_eqsgn_${v_ext}/hi-tsz+lens_final_v2_20mocks_eqsgn_${v_ext}_*/*.fits
    echo "Deleting ThumbStack filtarea, filtmask output..."
    trash /home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_final_v2_20mocks_eqsgn_${v_ext}/hi-tsz+lens_final_v2_20mocks_eqsgn_${v_ext}_*/output/thumbstack/*/tauring_filtarea.txt
    trash /home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_final_v2_20mocks_eqsgn_${v_ext}/hi-tsz+lens_final_v2_20mocks_eqsgn_${v_ext}_*/output/thumbstack/*/tauring_filtmask.txt
    trash-empty
    echo "Done cleaning and emptying trash."
done
