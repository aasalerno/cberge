function dz1 = modelP3(t,z,P)
% [original A. Salerno; re-hacked by C. Bergevin]
% defines ODEs to be solved by vdModelP#.m

global nR
persistent per per0 

% bookeeping re ode4
if P.solver==0
    % -----
    dz1 = zeros(P.N,1); % a column vector - for the solution
    ti = floor(t./P.stepSize);
    % -----
    if ti >= numel(P.L)
        ti = numel(P.L)-1;
    end
end

% -----
for j = 1:P.N
    
    % book keeping re the non-autonomous drive term (clunky as is)
    % (L1 acts on papilla, L2 on bundles)
    if P.solver==0
        if P.extD==0
            L1 = P.L(ti+1); L2= 0;
        elseif P.extD==1
            L1 = P.L(ti+1); L2= L1;
        end
    elseif P.solver==1
        if P.extT==0
            L1= 0; L2= L1;
        elseif P.extT==1
            L1= P.A*sin(2*pi*P.fe*t); 
            if (P.extD==0), L2=0;   else    L2= L1; end
        elseif P.extT==2
            L1= P.A*randn(1,1);
            if (P.extD==0), L2=0;   else    L2= P.A*randn(1,1); end
        end
    end
    
    % first oscillator will be the papilla, all others the 'bundles';
    % conditional for j>=2 is due to "nearest-neighbor" coupling at the ends
    if j == 1
        dz1(j) = z(j)*(i*P.wP+P.eP) +P.alpha*sum(P.dRP*(z(j)-z(2:end))) +P.beta*i*sum(P.dIP*(z(j)-z(2:end))) +L1;
    elseif j == 2
        dz1(j) = z(j)*(i*P.w(j)+P.e(j)-P.B(j)*abs(z(j))^2) +P.kap(j)*(P.dR(j)+i*P.dI(j))*(z(j+1)-z(j))+ P.alpha*P.dRP*(z(j)-z(1))+ P.beta*i*P.dIP*(z(j)-z(1)) +L2;
    elseif j == P.N
        dz1(j) = z(j)*(i*P.w(j)+P.e(j)-P.B(j)*abs(z(j))^2) +P.kap(j)*(P.dR(j)+i*P.dI(j))*(z(j-1)-z(j))+ P.alpha*P.dRP*(z(j)-z(1))+ P.beta*i*P.dIP*(z(j)-z(1)) +L2;
    else
        dz1(j) = z(j)*(i*P.w(j)+P.e(j)-P.B(j)*abs(z(j))^2) +P.kap(j)*(P.dR(j)+i*P.dI(j))*(z(j+1)+z(j-1)-2*z(j))+ P.alpha*P.dRP*(z(j)-z(1))+ P.beta*i*P.dIP*(z(j)-z(1)) +L2;
    end
    
end


% display some updating info to screen (kludgy, but works for both ode4 and
% ode45)
per0 = per;
per = (nR-1)/(4*(length(P.L)-1))*100;
nR = nR+1;
if (mod(fix(per),10)==0) & (fix(per)>fix(per0))
    %fprintf([num2str(fix(per)),'%% done \n'])
    fprintf([num2str(100*t/P.lengthT,2),'%% done \n'])
end

% bookeeping re ode45 (wants output as a column vector; be careful of .)
if (P.solver==1), dz1= dz1.';   end



