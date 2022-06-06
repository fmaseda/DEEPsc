function NormSCD=NormalizeRNAseq(SCD,method)
switch method
    case 'linear'
        % normalize by gene (columns); cell expressing gene at highest
        % level=1, all others in [0,1]
        NormSCD=SCD./(max(SCD)+eps());  
    case 'log'
        % log transform
        NormSCD=log(SCD+1);
    case {'logscale', 'loglinear'}
        % first do log transform, then scale linearly
        NormSCD=log(SCD+1)./(max(log(SCD+1))+eps());
end