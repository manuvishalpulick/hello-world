function main
clc 
clear all
% close all
%L_flat=60;
L_flat=15;
animationSkip=20;
seN=20;
N=2^9;
%N=2^10;
deltaX=L_flat/N;
x = 0:deltaX:L_flat;%+(2*deltaX);
ctimestep=2.75;
deltaT = deltaX^ctimestep; 
tt = seN*deltaT;
makeAnimation_det_fft(animationSkip, x, tt, L_flat,deltaX,N)
makeAnimation_det(animationSkip, x, tt, L_flat,deltaX,N)
end

function makeAnimation_det(animationSkip, x, tt, L_flat,deltaX,N)

newFolder = 'realization1';
str1 = strcat('.\',newFolder, '\*.txt');
textFiles = dir([str1]);
a   = cellfun(@num2str, struct2cell(textFiles), 'UniformOutput', false);
Out = sortrows(a.',6);              %sorted according to time so as to read in increasing times

b = length(textFiles);
t = (0:1:b).*tt;                    % corresponds to each time step, including the initial time step                    
formatSpec = 't = %0.3e,';          % to show time steps
t_string2 =  (sprintf(formatSpec,t)); % sprintf formats the data in t according to formatspec and returns the values into t_string2
t_new = strsplit(t_string2,',');    % splits str at delimeter ','

%% read the height profiles for every time steps
for k = 1:1:b   % 1:s:(b-last)
    str2 = strcat('.\', newFolder, '\', Out{k,1});
    c(:,k)  = dlmread(str2);
end
c(:,all(c==0)) = [];
c = [c(:,(1:end))];
q = max(size(t));

%% to specify animation conditions

v = VideoWriter('filmEvolution');
v.FrameRate = 5;  % Default 30
v.Quality = 100;    % Default 75
open(v)
hfig = figure;
l=length(c(:,1));
% x=linspace(0,L_flat,l);
% N1=length(x)
for i = 1:animationSkip:q
    Y = c(3:end-2,i);
    area(x',Y)
    ylim([0 5])
    xlim([0 L_flat])
    xlabel('x [-]','Fontsize',16)
    ylabel('h [-]','Fontsize',16)
%    ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
    legend(t_new(i))
    M(i) = getframe(hfig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards
end
M = M(~cellfun(@isempty,{M.cdata}));
writeVideo(v,M)
close(v)


end

function makeAnimation_det_fft(animationSkip, x, tt, L_flat,deltaX,N)

newFolder = 'realization1';
str1 = strcat('.\',newFolder, '\*.txt');
textFiles = dir([str1]);
a   = cellfun(@num2str, struct2cell(textFiles), 'UniformOutput', false);
Out = sortrows(a.',6);              %sorted according to time so as to read in increasing times

b = length(textFiles);
t = (0:1:b).*tt;                    % corresponds to each time step, including the initial time step                    
formatSpec = 't = %0.3e,';          % to show time steps
t_string2 =  (sprintf(formatSpec,t)); % sprintf formats the data in t according to formatspec and returns the values into t_string2
t_new = strsplit(t_string2,',');    % splits str at delimeter ','

%% read the height profiles for every time steps
for k = 1:1:b   % 1:s:(b-last)
    str2 = strcat('.\', newFolder, '\', Out{k,1});
    c(:,k)  = dlmread(str2);
end
ho=ones(size(x))+0.001*sin(6.*(x-30).^2);  % Defining the intitial condition once more
c(:,all(c==0)) = [];                        
c1 = [ho' c((3:end-2),:)];          % concantenating the initial condition with heights at other times
q = max(size(t));

%% to specify animation conditions

vid = VideoWriter('Fourier_evolution');
vid.FrameRate = 5;  % Default 30
vid.Quality = 100;    % Default 75
open(vid)
hfig = figure;
l=length(c(:,1));
%x=linspace(0,L_flat,l);
N1=length(x)
deltaX=(L_flat/N)

for i = 1:animationSkip:q
    Y = c1(:,i)
    Y_diff=Y-mean(Y);
    %% trying to anlayse in fourier domain
    hk=fft(Y_diff);
%    area(x',Y)
    hk2=fftshift(hk);
    P2 = abs(hk2/(N+1));
    P1 = P2(1:N/2+1); %since the middle value would be N/2 +1
    P1(2:end-1) = 2*P1(2:end-1)   % found in net dont know why
%    k= 2*pi.*(1./([2*deltaX:2*deltaX:L_flat]))';
    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L_flat))' ;
    S=P1.^2;
%    semilogy((S),'-','linewidth',2);
    semilogy(f,(S),'-','linewidth',2);
 %   semilogy(f,(S),'-','linewidth',2);
    xlabel('k [-]','Fontsize',16)
    ylabel('(S_k)[-]','Fontsize',16)
 
 %   ylim([0 5])
 %   xlim([0 L_flat])
    xlabel('k [-]','Fontsize',16)
    ylabel('Sk [-]','Fontsize',16)
%    ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
    legend(t_new(i))
    F(i) = getframe(hfig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards
end
F = F(~cellfun(@isempty,{F.cdata}));
writeVideo(vid,F)
close(vid)


end