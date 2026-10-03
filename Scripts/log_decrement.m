% fbg_mode1: table 160x7 -> [TIME UX UY UZ ROTX ROTY ROTZ]

M = table2array(fbg_mode1_1);
phi_mode_fbg = table2array(phi_mode_fbg);
phi_mode_fbg = phi_mode_fbg';
phi_mode_fbg = reshape(phi_mode_fbg,[48,1]);

nSens  = 8;
nSteps = size(M,1)/nSens;

% separa i sensori: FBG{i} = nSteps x 7
FBG = cell(nSens,1);
for i = 1:nSens
    FBG{i} = M(i:nSens:end, :);
end

% costruisci U globale: (6*nSens) x nSteps
U = zeros(6*nSens, nSteps);
for i = 1:nSens
    U((i-1)*6 + (1:6), :) = FBG{i}(:,2:7).';   % UX..ROTZ
end

% phi globale dai dati iniziali u(0)
phi = phi_mode_fbg(:,1);                                  % (6*nSens) x 1

% q(t) globale (1 x nSteps)
q = (phi.' * U) / (phi.' * phi);


time = FBG{1}(:,1);                             % nSteps x 1
q    = q.';                                     % nSteps x 1
