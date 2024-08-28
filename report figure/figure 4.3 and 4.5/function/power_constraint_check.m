function power_constraint_check(P_c,P_k,SNR,Num_U)
   power_c = trace(P_c'*P_c);
   power_k = zeros(1,Num_U);
   for i = 1:Num_U
        power_k(i) = trace(P_k(:,:,i)'*P_k(:,:,i)); 
   end

   if( power_c + sum(power_k) - SNR > 1e-4)
       error("Power constraint problem, please check the precoder power")
   end
end