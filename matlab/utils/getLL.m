function [LL,x,v,flag] = getLL(theta, Yobs,cfg,joint_ll,use_minimal_system)
is_for_standard_errors=true;
 [negLL, ~, ~,ll] = ...
     estimate_theta2(Yobs,cfg,theta,use_minimal_system,is_for_standard_errors);
 if joint_ll
     LL = -negLL;
 else
     LL=ll;
 end
end