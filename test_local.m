% Test local.

local_wave = make_local_waves(500);
local_emd_compare(2000,local_wave,700:-100:300,@benp_emd_local)

local_wave = make_local_triangle_waves(500);
local_emd_compare(2000,local_wave,700:-100:300,@benp_emd_local)

local_wave = make_transition_wave(5000);
local_emd_compare(2000,wave,700:-100:300,@benp_emd_local)