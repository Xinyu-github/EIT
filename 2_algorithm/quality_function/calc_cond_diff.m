function cond_diff=calc_cond_diff(ref_grid,recons_grid,n_ref)
%ss66
cond_diff=sum(abs(ref_grid(~isnan(ref_grid))-recons_grid(~isnan(recons_grid))))/n_ref;

end