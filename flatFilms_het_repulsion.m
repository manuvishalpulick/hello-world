
function t_rupt = flatFilms_het_repulsion(L,N,c,Tmp,gx,h_adjusted,A,p,endTime,seN);


format long g


tic

%% equation
%                                                                            ___
%  d H      d  /   3 d^3 H     1 dH       B   dH    \       __    d  /    /           \
%  --- = - ---|   H  ----- +   ----- -  ------ --    |  + _/2T   ---|    / H^3 N(X,T)  |
%  d T     d X \     d X^3     H dX       H^2  dX   /           d X \ \/              /


%% boundary conditions

% periodic

%% discretization
%           j      j 
%       2. H^2  * H^2                      _            _   _           _
%          i-1     i                    |   1       1  | |   j    j    |
%h1_l = _------------ _     h2_l = 0.5* |  ---  +  --- | |  H  - H     | 
%      |   j      j    |                |   j       j  | |_  i    i-1 _|
%      |  H    + H     |                |  H       H   |
%      |_  i-1    i   _|                |_  i-1     i _|  

%                 _                _   _           _
%                |   h_m     h_m    | |   j    j    |
%  h3_l =B* 0.5* |  ---  +  ---     | |  H  - H     | 
%                |     j         j  | |_  i    i-1 _|
%                |  H^2       H^2   |
%                |_    i-1       i _|  

%                 _                _   _           _
%                |   1       1      | |   j    j    |
%  h3_r = B*0.5* |  ---  +  ---     | |  H  - H     | 
%                |     j         j  | |_  i+1  i   _|
%                |  H^2       H^2   |
%                |_    i        i+1_|  

%           j      j 
%       2. H^2  * H^2                      _            _   _           _
%          i+1     i                    |   1       1  | |   j      j  |
%h1_r = _------------ _     h2_r = 0.5* |  ---  +  --- | |  H    - H   |
%      |   j      j    |                |   j       j  | |_  i+1    i _|
%      |  H    + H     |                |  H       H   |
%      |_  i+1    i   _|                |_  i+1     i _|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now here is the final discretization
%  _              
% |         j+1                             j+1                              j+1  
% | p2*h2_l* H        -p2*(3*h1_l + h1_r)* H        ( 3*p2*(h1_1+h1_r) + 1 ) H     ... 
% |_         i-2                            i-1                               i
%                                          _
%                     j+1             j+1   |                
%   -p2*(3*h1 + h2)* H        p2*h1* H      |    =  LHS --> A matrix
%                     i+1             i+2  _|                       


% _                                                                               _          
%| 0                                                                               | 
%| 0                                                                               | 
%|                                                                                 | 
%|previous        disj pressure                   noise term                       |     
%| soln             vdW only                                                       |
%|  j                                                                              |
%| H       -     p1 (h2_r - h2_l)  +    p3*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2)
%|  i                                                                              |        
%|                                                                                 |  
%|0                                                                                |  
%|_0                                                                              _| 


deltaX = L/N;                 % grid size
x = 0:deltaX:L;               % domain of the simulation
deltaT = deltaX^c            % time step size
flag = 0;                     % counter for storing the data files
p1 = deltaT/(deltaX^2);       % p1 - used in the explicit part
p2 = deltaT/(deltaX^4);       % p2 - used in the implicit part
p3 = 1/(deltaX)*sqrt(2*deltaT*Tmp);    % p3 - used for the noise term: Important to realize that factor 3 is not present in the avg mobility, if you do the proper scaling, so no need to have it in the sqrt here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%----------- Setting up the initial conditions ------------%%%
%P=(N/6*deltaX);
%h = ones(size(x));                        %  flat initial film
%%
%h=ones(size(x))+0.001*sin(6.*(x-30).^2); % for L=60
h=ones(size(x))+0.001*sin(6.*(x-7).^2);  %for L=15
% introducing an initial perturbation whichexcites multiple frequencies

h = [h(end-1), h(end), h, h(1), h(2)]';   % add two ghost points each side, and ensure periodicity employing periodic boundary conditions
Y_diff=h(3:h_adjusted-2)- mean(h(3:h_adjusted-2));
    hk=fft(Y_diff);
    hk2 = abs(hk/(L)); % Taking only the real values of the fourier transform and normalizing
    %which can still give us the required info
    P2=fftshift(hk2);
    P1 = P2(1:(N/2)+1); %since the middle value would be N/2 +1
    P1(2:end-1) = 2*P1(2:end-1);   %Have to double the magnitude as we 
    %ignored one side of the domain in the previous step
    %why is it 2:end-1 and not 1:end --> found in net dont know why
%    k= ((N/L).*(linspace(0,(N)/2,(N)/2+1))/N)';
    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L))';
%    f =(2*pi.*(linspace(0,(N+1),(N+1))./L))'
    S=P1.^2;
    semilogy(f(2:end-1),S(2:end-1),'-','linewidth',2);
    xlabel('k [-]','Fontsize',16)
    ylabel('(S_k) [-]','Fontsize',16)
figure
    plot(x,h(3:h_adjusted-2),'*')
ylim([0.998 1.002]);
xlim([0 L])
xlabel('x')
ylabel('h')
title('initial height profile')
%figure
%%
%h = [h(end-1), h(end), h, h(1), h(2)]';   % add two ghost points each side, and ensure periodicity employing periodic boundary conditions
InitArea = trapz(x,h(3:end-2));           % initial area of the film (mass) Just for confirmation

%%  preallocation

h1_r = zeros(size(h));
h1_l = zeros(size(h));
h2_r = zeros(size(h));
h2_l = zeros(size(h));
h3_r = zeros(size(h));
h3_l = zeros(size(h));
R = zeros(size(x));
noi1 = zeros(size(x));
noi2 = zeros(size(x));
B=0.1;  % constant of Part of disjoining pressure
%B=0; % For No repulsion case
k = linspace(1,h_adjusted,h_adjusted);   % index used in vectorizing the code; see the solver
e=0.1;
m=2*pi/(204);            %heterogeneity with a period of 100 nodes--> 10 P_het in L 
% For L=80, we have around 10 dominant wavelengths on it so Phet around
% 1024/10--> 100
%% Solver

rng shuffle                      % change the seed for every realization
for t = deltaT:deltaT:endTime    % time marching 
%     clear b
    R = randn(N+1,1);            % random numbers for each grid point (changes every time)
    gx_f = gx.*R;                % use the gx matrix as per the model
    g = sum(gx_f);               % Q-Wiener process definition as outlined in Grun et al, Diez et al, Lord et al
    g = [g(end-1); g(end); g'; g(1); g(2)];         % periodicity also in noise
    noi1 = (g(3:h_adjusted-2) + g(4:h_adjusted-1))./2;      % for the main domain
    noi2 = (g(2:h_adjusted-3) + g(3:h_adjusted-2))./2;      % for the main domain
    
    %% generate the pentagonal band in the sparse matrix
    
    % schemes outlined in Diez et al (2000), Grun et al (2004)[in appendix]
    % to discretize h^3 term
    h1_r = 2*h(k(3:h_adjusted-2)).^2.*h(k(4:h_adjusted-1)).^2./(h(k(3:h_adjusted-2))+h(k(4:h_adjusted-1)));          % see the discretization

    h1_l = 2*h(k(3:h_adjusted-2)).^2.*h(k(2:h_adjusted-3)).^2./(h(k(3:h_adjusted-2))+h(k(2:h_adjusted-3)));          % see the discretization

    h2_r =  (0.5.*(1+e*cos(m*(k(3:h_adjusted-2)))))'.*(1./h(k(3:h_adjusted-2)) + 1./h(k(4:h_adjusted-1))).*(h(k(4:h_adjusted-1)) - h(k(3:h_adjusted-2))); % see the discretization

    h2_l =  (0.5.*(1+e*cos(m*(k(3:h_adjusted-2)))))'.*(1./h(k(3:h_adjusted-2)) + 1./h(k(2:h_adjusted-3))).*(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3))); % see the discretization

    h3_r =  -0.5.*(1+e*cos(m*(k(3:h_adjusted-2))))'.*B.*(1./(h(k(3:h_adjusted-2))).^2 + 1./(h(k(4:h_adjusted-1))).^2).*(h(k(4:h_adjusted-1)) - h(k(3:h_adjusted-2))); % see the discretization

    h3_l =  -0.5.*(1+e*cos(m*(k(3:h_adjusted-2))))'.*B.*(1./(h(k(3:h_adjusted-2))).^2 + 1./(h(k(2:h_adjusted-3))).^2).*(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3))); % see the discretization

    
%      h3_r=(1./h(k(3:h_adjusted-2)) + 1./h(k(2:h_adjusted-3)));
%      h4_r=(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3)));
%      h2_2=h3_r.*h4_r.*h3_r;

    A(p) =[h1_l*p2; -(h1_r+3*h1_l)*p2; 3*(h1_r+h1_l)*p2+1; -(3*h1_r+h1_l)*p2; h1_r*p2];   % see the discretization - Arranging the 5 element space in the middle 
    %of the matrix through vectorization- p was already defined, each element in the RHS is of size=N+1 and p is of size N+1*5 (5--> -3,-2,-1,0,1) 
    
    %% generate the b-vector in Ax = b
    
    b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l+ h3_r- h3_l) + p3.*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2); 0; 0];                  % see the discretization
    % First and last 2 elements will be 0 as its is for the periodic
    % boundary conditions (eg: h(end)-h(1)=0 and so on)
%    b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2); 0; 0]; 
    h = A\b;                                                                                                           % solve for h
    h_min=min(h(:))
     if min(h(:)) <= 0.0001    % if the film height goes below a certain height stop the realization; it is insensitive to a value below 0.1. Rapid thinning
        break; % here we break out if the simulation and call it the rupture event and the time at this moment is the rupture time
     end
%====================================================================== 
%               Save data after every seN time steps
%======================================================================
    flag = flag + 1;                                   % count the no of time steps till seN time steps
    if(flag==seN)
        flag=0;
        dlmwrite(sprintf('Data_%1.8f.txt',t),h,'precision','%.16f','delimiter','\t')   % save it as a txt file to be used in matlab post processing script
% \t refers to a tab as a delimiter (default delimiter is a comma) , precision is also specified for accuracy and to prevent large space 
% consumption for too much precision  \

%some check to ensure if any realization has given singular matrix      
%         if(rcond(full(A)) < 1e-12)
%             break
%         else
%             continue
%         end
    end %end for if statement
%     Y_diff=h(3:h_adjusted-2)- mean(h(3:h_adjusted-2));
%     hk=fft(Y_diff);
%     hk2 = abs(hk/(N)); % Taking only the real values of the fourier transform and normalizing
%     %which can still give us the required info
%     P2=fftshift(hk2);
%     P1 = P2(1:(N/2)+1); %since the middle value would be N/2 +1
%     P1(2:end-1) = 2*P1(2:end-1);   %Have to double the magnitude as we 
%     %ignored one side of the domain in the previous step
%     %why is it 2:end-1 and not 1:end --> found in net dont know why
% %    k= ((N/L).*(linspace(0,(N)/2,(N)/2+1))/N)';
% %    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L))';
% %    f =(2*pi.*(linspace(0,(N+1),(N+1))./L))'
%     S=P1.^2;
%     semilogy(f,(S),'-','linewidth',2);
%     xlabel('$k [-]$','Fontsize',16)
%     ylabel('$(S_k)\:[-]$','Fontsize',16)
end
t_rupt = t;         % rupture time obtained from this simulation

toc

end

