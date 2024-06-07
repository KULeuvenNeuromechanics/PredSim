
(revisar bé tot)
% synergies only mesh points
% add track synergy weights



Settings: S.subject.synergies 

S.Syn -> S.subject.synergies


% Veure si canvio alguns noms (perquè siguin més clars / més intuïtius) i també on es guarden
% Posar-ho tot dins de S.Subject.synergies? O no aquí dins, sinó en "paral·lel"? A veure)
S.NSyn_r
S.NSyn_l
W.Syn_constr (in W or in Synergies?)
W.knownSynW (in W or in Synergies?)
S.knownSynW_idx
S.SynConstrLower
S.SynConstrUpper



Main changes in the code:

New variables:


New terms in the cost function:



New constraints:




Type of simulations that you can run:



TO DO:
- crec que no calen els parèntesis dels if's
- Better naming of settings
- Initialise settings for synergies
- automatic reading / indexing
- revise all cases are possible (symm, no symm, etc)
- reconstruction! (% TO DO: add muscle synergies term. For now: muscle synergies term will be the difference)