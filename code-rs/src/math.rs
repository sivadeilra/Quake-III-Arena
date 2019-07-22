pub fn Square<N>(n: N) -> N
where
    N: Copy,
    N: std::ops::Mul<N, Output = N>,
{
    n * n
}

// q_math
