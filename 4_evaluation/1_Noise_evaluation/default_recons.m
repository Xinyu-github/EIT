ff = fmdl;
ff.fwd_model = fmdl;
ff.jacobian_bkgnd.value = 1;
ff.type = 'inv_model';
ff.RtR_prior = 'eidors_default';
ff.hyperparameter.value = 0.05;
imgS = inv_solve_diff_GN_one_step(ff,v1,v2)

elem = imgS.elem_data-min(imgS.elem_data);
elem = elem/max(elem);

elem = elem.*(max(abs(img_diff.elem_data)));
elem = elem+min(img_diff.elem_data);

imgS.elem_data = elem;


figure, show_phase_contour_slice(img_diff,64,50,cm);
figure, show_phase_contour_slice(imgS,64,50,cm);


