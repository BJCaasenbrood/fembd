function transform(b, R, t)
%TRANSFORM Rigidly move the body: x -> R*x + t.
    b.X = b.X*R' + t(:)';
    b.Ref = b.Ref*R' + t(:)';
    if ~isempty(b.MassProps)
        b.MassProps.c = b.MassProps.c*R' + t(:)';
        b.MassProps.J = R*b.MassProps.J*R';
    end
end
