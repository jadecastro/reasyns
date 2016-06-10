function [vec, dvec] = buildBinaryCube(vec, dvec, vecSize)

if length(vec) == vecSize
    dvec = [dvec; vec];
    return
else
    vec = [vec -1];
    for j = -1:2:1
        vec(end) = j;
        [vec, dvec] = buildBinaryCube(vec, dvec, vecSize);
    end
    vec = vec(1:end-1);
end
