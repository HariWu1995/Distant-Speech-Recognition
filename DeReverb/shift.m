function O = shift(I,t)
% function O = shift(I,t)

if t<0
    O = [I(:,-t+1:end) zeros(size(I,1),-t) ];
else
    O = [zeros(size(I,1),t) I(:,1:end-t) ];
end