# p5-cubic-bezier

Cubic bezier with p5

## Usage

```java
final int SAMPLE = 500;

void setup () {
  size(500, 500);
  background(0);
  stroke(-1);

  CubicBezier bezier = new CubicBezier(0.42, 0, 0.58, 1);

  for (int i = 0; i < SAMPLE; i++) {
    float x = (float) ((double) width / SAMPLE * i);
    float y = (float) (height - bezier.solve((double) i / SAMPLE) * height);
    point(x, y);
  }
}
```
