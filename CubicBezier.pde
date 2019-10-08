class CubicBezier {
  static final int CUBIC_BEZIER_SPLINE_SAMPLES = 11;
  static final int K_MAX_NEWTON_TERATIONS = 4;
  static final double K_BEZIER_EPSILON = 1e-7;

  private double _ax, _ay, _bx, _by, _cx, _cy;
  private double _startGradient, _endGradient;
  private double _rangeMin = 0, _rangeMax = 1;
  private double[] _splineSamples = new double[CUBIC_BEZIER_SPLINE_SAMPLES];
  
  public CubicBezier (int p1x, int p1y, int p2x, int p2y) {
    this((double) p1x, (double) p1y, (double) p2x, (double) p2y);
  }
  
  public CubicBezier (float p1x, float p1y, float p2x, float p2y) {
    this((double) p1x, (double) p1y, (double) p2x, (double) p2y);
  }

  public CubicBezier (double p1x, double p1y, double p2x, double p2y) {
    this.initCoefficients(p1x, p1y, p2x, p2y);
    this.initGradients(p1x, p1y, p2x, p2y);
    this.initRange(p1y, p2y);
    this.initSpline();
  }

  private double sampleCurveX (double t) {
    return ((this._ax * t + this._bx) * t + this._cx) * t;
  }

  private double sampleCurveY (double t) {
    return ((this._ay * t + this._by) * t + this._cy) * t;
  }

  private double sampleCurveDerivativeX (double t) {
    return (3.0 * this._ax * t + 2.0 * this._bx) * t + this._cx;
  }

  private double sampleCurveDerivativeY (double t) {
    return (3.0 * this._ay * t + 2.0 * this._by) * t + this._cy;
  }

  public double getDefaultEpsilon () {
    return K_BEZIER_EPSILON;
  }

  private double solveCurveX (double x, double epsilon) {
    double t0 = 0, t1 = 0, t2 = x, x2 = 0, d2 = 0;
    int i;
    double deltaT = 1.0 / (CUBIC_BEZIER_SPLINE_SAMPLES - 1);
    for (i = 1; i < CUBIC_BEZIER_SPLINE_SAMPLES; i++) {
      if (x <= this._splineSamples[i]) {
        t1 = deltaT * i;
        t0 = t1 - deltaT;
        t2 = t0 + (t1 - t0) * (x - this._splineSamples[i - 1]) / (this._splineSamples[i] - this._splineSamples[i - 1]);
        break;
      }
    }

    double newtonEpsilon = Math.min(K_BEZIER_EPSILON, epsilon);
    for (i = 0; i < K_MAX_NEWTON_TERATIONS; i++) {
      x2 = this.sampleCurveX(t2) - x;
      if (Math.abs(x2) < newtonEpsilon) return t2;
      d2 = this.sampleCurveDerivativeX(t2);
      if (Math.abs(d2) < K_BEZIER_EPSILON) break;
      t2 = t2 - x2 / d2;
    }
    if (Math.abs(x2) < epsilon) return t2;

    while (t0 < t1) {
      x2 = this.sampleCurveX(t2);
      if (Math.abs(x2 - x) < epsilon) return t2;
      if (x > x2) t0 = t2; 
      else t1 = t2;
      t2 = (t0 + t1) * 0.5;
    }
    return t2;
  }

  public double solve (double x) {
    return this.solveWithEpsilon(x, K_BEZIER_EPSILON);
  }

  public double solveWithEpsilon (double x, double epsilon) {
    if (x < 0.0) return 0.0 + this._startGradient * x;
    if (x > 1.0) return 1.0 + this._endGradient * (x - 1.0);
    return this.sampleCurveY(this.solveCurveX(x, epsilon));
  }

  public double slope (double x) {
    return this.slopeWithEpsilon(x, K_BEZIER_EPSILON);
  }

  public double slopeWithEpsilon(double x, double epsilon) {
    x = Math.max(0.0, Math.min(x, 1.0));
    double t = this.solveCurveX(x, epsilon);
    double dx = this.sampleCurveDerivativeX(t);
    double dy = this.sampleCurveDerivativeY(t);
    return dy / dx;
  }

  public double getX1 () {
    return this._cx / 3.0;
  }

  public double getY1 () {
    return this._cy / 3.0;
  }

  public double getX2 () {
    return (this._bx + this._cx) / 3.0 + this.getX1();
  }

  public double getY2 () {
    return (this._by + this._cy) / 3.0 + this.getY1();
  }

  public double getRangeMin() {
    return this._rangeMin;
  }

  public double getRangeMax() {
    return this._rangeMax;
  }

  private void initCoefficients (double p1x, double p1y, double p2x, double p2y) {
    this._cx = 3.0 * p1x;
    this._cy = 3.0 * p1y;
    this._bx = 3.0 * (p2x - p1x) - this._cx;
    this._by = 3.0 * (p2y - p1y) - this._cy;
    this._ax = 1.0 - this._cx - this._bx;
    this._ay = 1.0 - this._cy - this._by;
  }

  private void initGradients (double p1x, double p1y, double p2x, double p2y) {
    if (p1x > 0) this._startGradient = p1y / p1x;
    else if (p1y == 0 && p2x > 0) this._startGradient = p2y / p2x;
    else if (p1y == 0 && p2y == 0) this._startGradient = 1;
    else this._startGradient = 0;

    if (p2x < 1) this._endGradient = (p2y - 1) / (p2x - 1);
    else if (p2y == 1 && p1x < 1) this._endGradient = (p1y - 1) / (p1x - 1);
    else if (p2y == 1 && p1y == 1) this._endGradient = 1;
    else this._endGradient = 0;
  }

  private void initRange (double p1y, double p2y) {
    if (0 <= p1y && p1y < 1 && 0 <= p2y && p2y <= 1) return;

    double a = this._ay * 3.0;
    double b = this._by * 2.0;
    double c = this._cy;
    if (Math.abs(a) < K_BEZIER_EPSILON && Math.abs(b) < K_BEZIER_EPSILON) return;

    double t1 = 0, t2 = 0;
    if (Math.abs(a) < K_BEZIER_EPSILON) {
      t1 = -c / b;
    } else {
      double discriminant = b * b - 4 * a * c;
      if (discriminant < 0) return;
      double discriminantSqrt = Math.sqrt(discriminant);
      t1 = (-b + discriminantSqrt) / (2 * a);
      t1 = (-b - discriminantSqrt) / (2 * a);
    }

    double sol1 = 0, sol2 = 0;
    if (0 < t1 && t1 < 1) sol1 = this.sampleCurveY(t1);
    if (0 < t2 && t2 < 1) sol2 = this.sampleCurveY(t2);

    this._rangeMin = Math.min(this._rangeMin, Math.min(sol1, sol2));
    this._rangeMax = Math.max(this._rangeMax, Math.max(sol1, sol2));
  }

  private void initSpline () {
    double deltaT = 1.0 / (CUBIC_BEZIER_SPLINE_SAMPLES - 1);
    for (int i = 0; i < CUBIC_BEZIER_SPLINE_SAMPLES; i++) {
      this._splineSamples[i] = this.sampleCurveX(i * deltaT);
    }
  }
}
