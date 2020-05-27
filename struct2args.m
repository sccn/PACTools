% Transform strcuture to  pairs 'argname', 'argvalue'
function opt = struct2args(g)
fields = fieldnames(g);
values = struct2cell(g);
params = { fields{:}; values{:} };
options = [ vararg2str( { params{:} } ) ];
opt     = eval(['{' options '}']);