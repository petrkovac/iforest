# Isolation Forest
Implementation based on: http://cs.nju.edu.cn/zhouzh/zhouzh.files/publication/icdm08b.pdf
## Installation
Just include `isolation_forest.h` in your project and you are ready to go.
## Requirements
Requires a C++11 or newer compiler:
- Windows: MSVC++ 2013 or newer
- Linux: gcc 4.9.1+ (compile with `--std=c++11`)

## Features
- easy integration, easy to use
- low memory requirements
- fast
- supports serialization/deserialization

## Example
```c++
#include <iostream>
#include <stdint.h>
#include <random>
#include <vector>
#include <string>
#include "isolation_forest.h"

int main(int argc, char* argv[])
{
	std::mt19937 rng(12345);
	std::vector<std::array<uint32_t, 4>> data;

	// generate some random 4D datapoints
	for (uint32_t i = 0; i < 1000; i++)
	{
		std::array<uint32_t, 4> temp;
		for (uint32_t j = 0; j < 4; j++)
		{
			temp[j] = iforest::UniformRandomNumber<uint32_t>::GenerateNext(rng, 0, 1000);
		}
		data.push_back(temp);
	}

	// add a few anomalies
	for (uint32_t i = 0; i < 5; i++)
	{
		std::array<uint32_t, 4> temp;
		for (uint32_t j = 0; j < 4; j++)
		{
			temp[j] = iforest::UniformRandomNumber<uint32_t>::GenerateNext(rng, 800, 1500);
		}
		data.push_back(temp);
	}

	iforest::IsolationForest<uint32_t, 4> forest;

	if (!forest.Build(50, 12345, data, 100))
	{
		std::cerr << "Failed to build Isolation Forest.\n";
		return 1;
	}

	std::vector<double> anomaly_scores;

	if (!forest.GetAnomalyScores(data, anomaly_scores))
	{
		std::cerr << "Failed to calculate anomaly scores.\n";
		return 2;
	}

	for (uint32_t i = 995; i < 1005; i++)
	{
		std::cout << "Anomaly_score[" << i << "] " << anomaly_scores[i] << "\n";
	}

	return 0;
}
```

### Example output
```
Anomaly_score[995] 0.487166
Anomaly_score[996] 0.457642
Anomaly_score[997] 0.497074
Anomaly_score[998] 0.475812
Anomaly_score[999] 0.468102
Anomaly_score[1000] 0.67561
Anomaly_score[1001] 0.665234
Anomaly_score[1002] 0.695242
Anomaly_score[1003] 0.6442
Anomaly_score[1004] 0.678004
```

## License
The class is licensed under the [MIT License](http://opensource.org/licenses/MIT):

Copyright © 2016 Peter Kovac

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.