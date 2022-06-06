function Noisyobj=AddNoise(obj,intensity,lowerlim,upperlim)
	Noisyobj=obj+intensity*randn(size(obj));
	Noisyobj=max(Noisyobj,lowerlim);
	Noisyobj=min(Noisyobj,upperlim);
end
