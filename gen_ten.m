function X = gen_ten(U)

s = cellfun(@(x)size(x,1),U).';
U = {U{end:-1:1}};
X = reshape(sum(khatrirao_Z(U),2),s);

