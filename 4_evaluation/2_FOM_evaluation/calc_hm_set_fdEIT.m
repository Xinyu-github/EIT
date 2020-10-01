function hm_img = calc_hm_set_fdEIT(img,frac)
% hm_img= CALC_HA_SET(img)
% hm_img is szxszxNimg. It is 1 inside the Half Ampl Set
% frac is the fraction of maximum (0.5 or 0.25)
% hm_img expects conductive changes. Use calc_hm_set(-img,frac) for non-c

% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id: calc_hm_set.m 4809 2015-03-29 11:55:10Z bgrychtol-ipa $

sz = size(img);
[x,y]=meshgrid(linspace(-1,1,sz(1)),linspace(-1,1,sz(2))); map = x.^2+y.^2<1.1;

hm_img = logical(zeros(size(img)));
for i=1:size(img,3);
   imi = img(:,:,i); imi= imi(map);

   hmi= logical(zeros(sz));
   hmi(map) = imi <= (mean(imi(imi~=0)) * frac) & imi ~= 0;
   hm_img(:,:,i) = hmi;
end
