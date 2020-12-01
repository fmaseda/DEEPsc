function [accuracy,precision]=CorrErrorKnownResult(Corr)
    % accuracy measures whether or not Corr=1 when it should (i.e. on the
    % diagonal)
	accuracy=mean(1-diag(Corr));
    % precision measures whether or not Corr=1 ONLY when it should (i.e.
    % the sum over positions is 1)
    precision=mean(abs((1-sum(Corr,2))./(size(Corr,2)-1)));
end