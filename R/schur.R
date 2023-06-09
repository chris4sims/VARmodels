#' Schur decomposition
#'
#' A front end for lapack complex schur decomposition `zgees`.
#' @param A The matrix to be decomposed
#'
#' @export
#' @md
#' 
schur <- function(A) {
    ## Matrix::Schur would provide a more platform-independent approach.
    ## But to use schdiv(), would need to post-process Schur output to
    ## get full complex, fully triangular, T.
  if ( !is.loaded("zgees")) dyn.load(La_library, local=TRUE)
  if (is.null(dim(A))) dim(A) <- c(1,1)
  if ( dim(A)[1] ==0) return(list(X=matrix(0,0,0), info=0))
  n <- dim(A)[1]
  n <- as.integer(n)
  sdim <- as.integer(0)
    schout <- .Fortran("zgees", "V","N","dum", n, T=as.complex(A), n, sdim,
                       ev=complex(n), Q=complex(n*n), n,
                       work=complex( 33 * n ), LWORK=as.integer(33 * n), rwork = double(n),
                       bwork = "dum", info=integer(1))
  return(list(Q=matrix(schout$Q,n), T=matrix(schout$T,n)))
  ##       SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
  ##      $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
  ## *
  ## *  -- LAPACK driver routine (version 3.1) --
  ## *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  ## *     November 2006
  ## *
  ## *     .. Scalar Arguments ..
  ##       CHARACTER          JOBVS, SORT
  ##       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
  ## *     ..
  ## *     .. Array Arguments ..
  ##       LOGICAL            BWORK( * )
  ##       DOUBLE PRECISION   RWORK( * )
  ##       COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
  ## *     ..
  ## *     .. Function Arguments ..
  ##       LOGICAL            SELECT
  ##       EXTERNAL           SELECT
  ## *     ..
  ## *
  ## *  Purpose
  ## *  =======
  ## *
  ## *  ZGEES computes for an N-by-N complex nonsymmetric matrix A, the
  ## *  eigenvalues, the Schur form T, and, optionally, the matrix of Schur
  ## *  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
  ## *
  ## *  Optionally, it also orders the eigenvalues on the diagonal of the
  ## *  Schur form so that selected eigenvalues are at the top left.
  ## *  The leading columns of Z then form an orthonormal basis for the
  ## *  invariant subspace corresponding to the selected eigenvalues.
  ## *
  ## *  A complex matrix is in Schur form if it is upper triangular.
  ## *
  ## *  Arguments
  ## *  =========
  ## *
  ## *  JOBVS   (input) CHARACTER*1
  ## *          = 'N': Schur vectors are not computed;
  ## *          = 'V': Schur vectors are computed.
  ## *
  ## *  SORT    (input) CHARACTER*1
  ## *          Specifies whether or not to order the eigenvalues on the
  ## *          diagonal of the Schur form.
  ## *          = 'N': Eigenvalues are not ordered:
  ## *          = 'S': Eigenvalues are ordered (see SELECT).
  ## *
  ## *  SELECT  (external procedure) LOGICAL FUNCTION of one COMPLEX*16 argument
  ## *          SELECT must be declared EXTERNAL in the calling subroutine.
  ## *          If SORT = 'S', SELECT is used to select eigenvalues to order
  ## *          to the top left of the Schur form.
  ## *          IF SORT = 'N', SELECT is not referenced.
  ## *          The eigenvalue W(j) is selected if SELECT(W(j)) is true.
  ## *
  ## *  N       (input) INTEGER
  ## *          The order of the matrix A. N >= 0.
  ## *
  ## *  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
  ## *          On entry, the N-by-N matrix A.
  ## *          On exit, A has been overwritten by its Schur form T.
  ## *
  ## *  LDA     (input) INTEGER
  ## *          The leading dimension of the array A.  LDA >= max(1,N).
  ## *
  ## *  SDIM    (output) INTEGER
  ## *          If SORT = 'N', SDIM = 0.
  ## *          If SORT = 'S', SDIM = number of eigenvalues for which
  ## *                         SELECT is true.
  ## *
  ## *  W       (output) COMPLEX*16 array, dimension (N)
  ## *          W contains the computed eigenvalues, in the same order that
  ## *          they appear on the diagonal of the output Schur form T.
  ## *
  ## *  VS      (output) COMPLEX*16 array, dimension (LDVS,N)
  ## *          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
  ## *          vectors.
  ## *          If JOBVS = 'N', VS is not referenced.
  ## *
  ## *  LDVS    (input) INTEGER
  ## *          The leading dimension of the array VS.  LDVS >= 1; if
  ## *          JOBVS = 'V', LDVS >= N.
  ## *
  ## *  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
  ## *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  ## *
  ## *  LWORK   (input) INTEGER
  ## *          The dimension of the array WORK.  LWORK >= max(1,2*N).
  ## *          For good performance, LWORK must generally be larger.
  ## *
  ## *          If LWORK = -1, then a workspace query is assumed; the routine
  ## *          only calculates the optimal size of the WORK array, returns
  ## *          this value as the first entry of the WORK array, and no error
  ## *          message related to LWORK is issued by XERBLA.
  ## *
  ## *  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
  ## *
  ## *  BWORK   (workspace) LOGICAL array, dimension (N)
  ## *          Not referenced if SORT = 'N'.
  ## *
  ## *  INFO    (output) INTEGER
  ## *          = 0: successful exit
  ## *          < 0: if INFO = -i, the i-th argument had an illegal value.
  ## *          > 0: if INFO = i, and i is
  ## *               <= N:  the QR algorithm failed to compute all the
  ## *                      eigenvalues; elements 1:ILO-1 and i+1:N of W
  ## *                      contain those eigenvalues which have converged;
  ## *                      if JOBVS = 'V', VS contains the matrix which
  ## *                      reduces A to its partially converged Schur form.
  ## *               = N+1: the eigenvalues could not be reordered because
  ## *                      some eigenvalues were too close to separate (the
  ## *                      problem is very ill-conditioned);
  ## *               = N+2: after reordering, roundoff changed values of
  ## *                      some complex eigenvalues so that leading
  ## *                      eigenvalues in the Schur form no longer satisfy
  ## *                      SELECT = .TRUE..  This could also be caused by
  ## *                      underflow due to scaling.
  ## *
  ## *  =====================================================================
}
