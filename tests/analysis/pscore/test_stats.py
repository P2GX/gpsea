import typing
import pytest

import numpy as np

from gpsea.analysis.pscore.stats import MannWhitneyStatistic, TTestStatistic


class TestMannWhitneyStatistic:
    @pytest.fixture(scope="class")
    def statistic(self) -> MannWhitneyStatistic:
        return MannWhitneyStatistic()

    @pytest.mark.parametrize(
        "x, y, expected",
        [
            (
                (
                    1.0,
                    2.0,
                    3.0,
                ),
                (
                    1.0,
                    2.0,
                    3.0,
                ),
                1.0,
            ),
            (
                (
                    11.0,
                    15,
                    8.0,
                    12.0,
                ),
                (
                    4.0,
                    2.0,
                    3.0,
                    3.5,
                    4.0,
                ),
                0.01945103333136247,
            ),
        ],
    )
    def test_compute_pval(
        self,
        statistic: MannWhitneyStatistic,
        x: typing.Sequence[float],
        y: typing.Sequence[float],
        expected: float,
    ):
        actual = statistic.compute_pval((x, y))

        assert actual.pval == pytest.approx(expected)

    def test_compute_pval__with_nan(
        self,
        statistic: MannWhitneyStatistic,
    ):
        x = (1.0, 2.0, 3.0, np.nan)
        y = (1.0, 2.0, 3.0, float("nan"))

        actual = statistic.compute_pval((x, y))

        assert actual.pval == pytest.approx(1.0)


class TestTTestStatistic:
    @pytest.fixture(scope="class")
    def statistic(self) -> TTestStatistic:
        return TTestStatistic()

    @pytest.mark.parametrize(
        "x, y, expected",
        [
            (
                (
                    1.0,
                    2.0,
                    3.0,
                ),
                (
                    1.0,
                    2.0,
                    3.0,
                ),
                1.0,
            ),
            (
                (
                    11.0,
                    15,
                    8.0,
                    12.0,
                ),
                (
                    4.0,
                    2.0,
                    3.0,
                    3.5,
                    4.0,
                ),
                0.0004749950471148506,
            ),
        ],
    )
    def test_compute_pval(
        self,
        statistic: TTestStatistic,
        x: typing.Sequence[float],
        y: typing.Sequence[float],
        expected: float,
    ):
        result = statistic.compute_pval((x, y))

        assert result.pval == pytest.approx(expected)

    def test_compute_pval__with_nan(
        self,
        statistic: TTestStatistic,
    ):
        x = (1.0, 2.0, 3.0, np.nan, np.nan)
        y = (1.0, 2.0, 3.0, float("nan"))

        result = statistic.compute_pval((x, y))

        assert result.pval == pytest.approx(1.0)
