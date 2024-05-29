from mango.mango_features.excluded_regions import RectangularPrism, Sphere
import numpy as np

if __name__ == '__main__':
    c1, c2 = np.array([0, 0, 0]), np.array([10, 10, 10])
    region = RectangularPrism(c1=c1, c2=c2)

    # Create sphere:
    region2 = Sphere(diameter=15, center=np.array([0, 0, 0]))

    print(region)
    print(region2)