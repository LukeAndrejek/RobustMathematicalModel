clear all

% This script implements the Arabidopsis leaf model
% This model involves seven variables: H, S, T, M, K, R, W
% H: HD-ZIP3, S: AS2, T: ta-siRNA, M: miR166, K: KANADI, R: ARF3, W: WOX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control how frequently the functions are plotted
plotGap = 58;

% Set noise parameters
noiseHC = 0;
noiseSC = 0;
noiseTC = 0;
noiseMC = 0;
noiseRC = 0;
noiseKC = 0;
noiseWC = 0;

% Set the standard deviation for noise truncation
noiseSD = 3;

% Multiply source and degradation parameters by a constant c
c = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter domain ranges
endT = 96;
endX = 20;
t1 = 72;

% Set the maximum value of dt
dtDefault = 0.01;

% Set the length of the cell wall for visual purposes
wallLength = 0.2;
dxWall = endX/6;

% Choose a constant near 0 to avoid division by 0 in some calculations
err = 1e-7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model Parameter Values

% Diffusion constants
DT = 9;
DM = 9;

% Degredation constants
dH = 0.15 * c;
dS = 0.15 * c;
dT = 0.1 * c;
dM = 0.1 * c;
dR = 0.15 * c;
dK = 0.15 * c;
dW = 0.15 * c;

% Source constants
sigmaH = 43 * dH * c;
sigmaS = 6 * dS * c;
sigmaT = 1.58;
sigmaM = 45.3;
sigmaR = 10 * dR * c;
sigmaK = 14 * dK * c;
sigmaW = 5 * dW * c;

% Interaction constants
kHM = 1;
kMH = kHM;
kTR = 1;
kRT = kTR;
kKW = 2;
kSM = 2;
kSR = 1;
kWS = 2;
kKRS = 2;

% Hill function parameters
hillExpM = 6;
hillExpT = 6;

% Set initial conditions
N = ones(6,1);
Z = zeros(6,1);
H0 = Z;
S0 = Z;
T0 = Z;
M0 = Z;
R0 = Z;
K0 = Z;
W0 = Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional setup

% Control how frequently the value of t is displayed
timerGap = 1;
plotBurst = 1.5*dtDefault;
timerBurst = 1.5*dtDefault;

% Define various spatial discretizations
dx = endX/6;
Xplot = [0; wallLength/2-err; wallLength/2+err];
for i = 1:5
    Xplot = [Xplot; i*dx-wallLength/2-err; i*dx-wallLength/2+err; ...
        i*dx+wallLength/2-err; i*dx+wallLength/2+err];
end
Xplot = [Xplot; endX-wallLength/2-err; endX-wallLength/2+err; endX];
XplotRNA = [0; wallLength/2];
for i = 1:5
    XplotRNA = [XplotRNA; i*dx-wallLength/2; i*dx+wallLength/2];
end
XplotRNA = [XplotRNA; endX-wallLength/2; endX];

% This vector helps with designing the protein diffusion operator
UiProtein = [1; 2; 2; 2; 2; 1];

% These vectors help plot functions
expand = [kron(eye(6), [0; 0; 1; 1]); zeros(2,6)];
expandRNA = [[1 0 0 0 0 0]; kron(eye(6),[1; 1]); [0 0 0 0 0 1]];

% Define supports for specific proteins
supportS = [1; 0; 0; 0; 0; 1];
supportM = [0; 0; 0; 0; 0; 1];
supportK = [0; 0; 0; 0; 1; 1];

% This section is used to make a truncated normal distribution
normDist = makedist('Normal');
truncDist = truncate(normDist,-noiseSD,noiseSD);
icdfN = 10000;
X2 = 0:1/(icdfN-1):1;
icdfApprox = icdf(truncDist, 0:1/(icdfN-1):1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main loop

t = 0;
while t < endT
    
    % Determine noise vectors
    noiseH = icdfApprox(ceil(icdfN*rand(6,1)))';
    noiseS = icdfApprox(ceil(icdfN*rand(6,1)))';
    noiseT = icdfApprox(ceil(icdfN*rand(6,1)))';
    noiseM = icdfApprox(ceil(icdfN*rand(6,1)))';
    noiseR = icdfApprox(ceil(icdfN*rand(6,1)))';
    noiseK = icdfApprox(ceil(icdfN*rand(6,1)))';
    noiseW = icdfApprox(ceil(icdfN*rand(6,1)))';
    
    % Create Hill functions
    if min(H0) == 0 || min(M0) == 0
        hillHM = (M0).^hillExpM ./ ...
            (H0.^hillExpM + (M0).^hillExpM + err*N);
    else
        hillHM = N ./ (((H0./M0)).^hillExpM + N);
    end
    
    if min(T0) == 0 || min(R0) == 0
        hillTR = (T0).^hillExpT ./ ...
            (R0.^hillExpT + (T0).^hillExpT + err*N);
    else
        hillTR = N ./ (((R0./T0)).^hillExpT + N);
    end
    
    % Adjust t to align with 
    dt = dtDefault;
    if t < t1
        dt = min(dt, t1 - t);
    else
        dt = min(dt, endT-t);
    end
    t = t + dt;
    
    % Define source functions
    ventral = M0 > H0;
    dorsal = N - ventral;
    supportM = [0; 0; 0; 0; 0; 1] .* (M0 >= H0);
    sourceH = sigmaH;
    sourceS = supportS .* (sigmaS ./ (N + kWS*W0 + kKRS*K0.*R0));
    sourceT = sigmaT * dorsal;
    sourceM = (sigmaM ./ (N + kSM*S0)) .* supportM;
    sourceR = sigmaR ./ (N + kSR*S0);
    sourceK = sigmaK * ventral .* supportK;
    sourceW = (t > t1) * (sigmaW ./ (N + kKW*K0)) .* ventral;
    
    % Apply operators
    H = H0 + dt * sourceH - dt*dH*H0 - dt*kHM*hillHM.*H0.*M0 ...
        + noiseHC*sqrt(dt)*noiseH.*H0;
    S = S0 + dt*sourceS - dt*dS*S0 + noiseSC*sqrt(dt)*noiseS.*S0;
    T = T0 + dt*sourceT - dt*dT*T0 + dt * dxWall^-2 * DT ...
        * ([T0(2:6); 0] - UiProtein.*T0 + [0; T0(1:5)]) ...
        - dt*kTR*hillTR.*T0.*R0 + noiseTC*sqrt(dt)*noiseT.*T0;
    M = M0 + dt*sourceM - dt*dM*M0 + dt * dxWall^-2 * DM ...
        * ([M0(2:6); 0] - UiProtein.*M0 + [0; M0(1:5)]) ...
        - dt*kMH*hillHM.*H0.*M0 + noiseMC*sqrt(dt)*noiseM.*M0;
    R = R0 + dt*sourceR - dt*dR*R0 - dt*kRT*hillTR.*T0.*R0 ...
        + noiseRC*sqrt(dt)*noiseR.*R0;
    K = K0 + dt*sourceK - dt*dK*K0 + noiseKC*sqrt(dt)*noiseK.*K0;
    W = W0 + dt*sourceW - dt*dW*W0 + noiseWC*sqrt(dt)*noiseW.*W0;
    
    % Stop the simulation if any functions are negative
    % This can occur if dt is too small
    if min([min(H), min(S), min(T), min(M), min(R), min(K), min(W)]) < 0
        disp('Negative function values')
        return
    end

    % Plot functions occasionally
    if mod(t,plotGap) < plotBurst && t > plotBurst

        clf
        subplot(1,3,1)
        plot(Xplot,expand*H, 'color', [1,0,0], 'LineWidth', 2);
        hold on
        plot(XplotRNA,expandRNA*M, 'color', [0,1,0], 'LineWidth', 2);
        hold off
        ylim([0 1.3*max([max(H), max(M)])]);
        xlim([0 endX]);
        lHM = legend('HD-ZIPIII', 'miR166', 'Location', 'northwest' );
        lHM.FontSize = 10;
        lHM.NumColumns = 2;
        set(gca,'xtick',endX/12:endX/6:endX-endX/12)
        set(gca,'xticklabel',{'{\it C}_1', '{\it C}_2', '{\it C}_3', ...
            '{\it C}_4', '{\it C}_5', '{\it C}_6',})
        xlabel('Cells')
        ylabel('TPM')
        set(gca,'FontSize',20)
        
        subplot(1,3,2)
        plot(XplotRNA,expandRNA*T, 'color', [1,0.5,0], 'LineWidth', 2);
        hold on
        plot(Xplot,expand*R, 'color', [0,0,1], 'LineWidth', 2);
        hold off
        ylim([0 1.3*max([max(T), max(R)])]);
        xlim([0 endX]);
        lTR = legend('tasiARF', 'ARF', 'Location', 'northwest' );
        lTR.FontSize = 10;
        lTR.NumColumns = 2;
        set(gca,'xtick',endX/12:endX/6:endX-endX/12)
        set(gca,'xticklabel',{'{\it C}_1', '{\it C}_2', '{\it C}_3', ...
            '{\it C}_4', '{\it C}_5', '{\it C}_6',})
        xlabel('Cells')
        ylabel('TPM')
        set(gca,'FontSize',20)
        clock = title(['TPM at t = ' num2str(round(t))]);
        clock.FontSize = 12;
        
        subplot(1,3,3)
        plot(Xplot,expand*S, 'color', [.9,.9,0], 'LineWidth', 2);
        hold on
        plot(Xplot,expand*K, 'color', [1,0,1], 'LineWidth', 2);
        plot(Xplot,expand*W, 'color', [0, 1, 1], 'LineWidth', 2);
        hold off
        ylim([0 max([max(S), max(K), max(W)])]);
        xlim([0 endX]);
        lSK = legend('AS', 'KANADI', 'WOX', 'Location', 'northwest' );
        lSK.FontSize = 10;
        set(gca,'xtick',endX/12:endX/6:endX-endX/12)
        set(gca,'xticklabel',{'{\it C}_1', '{\it C}_2', '{\it C}_3', ...
            '{\it C}_4', '{\it C}_5', '{\it C}_6',})
        xlabel('Cells')
        ylabel('TPM')
        set(gca,'FontSize',20)
        
        h=gcf;
        set(h,'PaperPositionMode','auto');
        set(h,'PaperOrientation','landscape');
        set(h,'Position',[50 50 1200 300]);
        
        % Pause briefly so that the figure can display
        pause(0.001)
        
    end
    
    % Display t occasionally
    if mod(t,timerGap) < timerBurst && t > timerBurst
        disp(['t: ' num2str(t)])
    end
    
    % Prepare for next loop
    H0 = H;
    S0 = S;
    T0 = T;
    M0 = M;
    R0 = R;
    K0 = K;
    W0 = W;

end

% Plot functions once the model has finished running
clf
subplot(1,3,1)
plot(Xplot,expand*H, 'color', [1,0,0], 'LineWidth', 2);
hold on
plot(XplotRNA,expandRNA*M, 'color', [0,1,0], 'LineWidth', 2);
hold off
ylim([0 1.3*max([max(H), max(M)])]);
xlim([0 endX]);
lHM = legend('HD-ZIPIII', 'miR166', 'Location', 'northwest' );
lHM.FontSize = 10;
lHM.NumColumns = 2;
set(gca,'xtick',endX/12:endX/6:endX-endX/12)
set(gca,'xticklabel',{'{\it C}_1', '{\it C}_2', '{\it C}_3', ...
    '{\it C}_4', '{\it C}_5', '{\it C}_6',})
xlabel('Cells')
ylabel('TPM')
set(gca,'FontSize',20)

subplot(1,3,2)
plot(XplotRNA,expandRNA*T, 'color', [1,0.5,0], 'LineWidth', 2);
hold on
plot(Xplot,expand*R, 'color', [0,0,1], 'LineWidth', 2);
hold off
ylim([0 1.3*max([max(T), max(R)])]);
xlim([0 endX]);
lTR = legend('tasiARF', 'ARF', 'Location', 'northwest' );
lTR.FontSize = 10;
lTR.NumColumns = 2;
set(gca,'xtick',endX/12:endX/6:endX-endX/12)
set(gca,'xticklabel',{'{\it C}_1', '{\it C}_2', '{\it C}_3', ...
    '{\it C}_4', '{\it C}_5', '{\it C}_6',})
xlabel('Cells')
ylabel('TPM')
set(gca,'FontSize',20)
clock = title(['TPM at t = ' num2str(round(t))]);
clock.FontSize = 12;

subplot(1,3,3)
plot(Xplot,expand*S, 'color', [.9,.9,0], 'LineWidth', 2);
hold on
plot(Xplot,expand*K, 'color', [1,0,1], 'LineWidth', 2);
plot(Xplot,expand*W, 'color', [0, 1, 1], 'LineWidth', 2);
hold off
ylim([0 max([max(S), max(K), max(W)])]);
xlim([0 endX]);
lSK = legend('AS', 'KANADI', 'WOX', 'Location', 'northwest' );
lSK.FontSize = 10;
set(gca,'xtick',endX/12:endX/6:endX-endX/12)
set(gca,'xticklabel',{'{\it C}_1', '{\it C}_2', '{\it C}_3', ...
    '{\it C}_4', '{\it C}_5', '{\it C}_6',})
xlabel('Cells')
ylabel('TPM')
set(gca,'FontSize',20)

h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 300]);