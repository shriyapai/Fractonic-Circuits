function[U]=phaseAction

u1 = randU(1); u2= randU(1); u3=randU(1);

U=sparse(diag([u1, u2, u3]));

