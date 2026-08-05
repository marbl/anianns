import csv

from anianns.union_find import SatelliteDSU


def test_satellite_dsu_unions_distal_links_transitively():
    satellites = SatelliteDSU()
    first = satellites.add_satellite("chr1", 100, 200, name="A")
    second = satellites.add_satellite("chr1", 500, 700, name="B")
    third = satellites.add_satellite("chr1", 1000, 1200, name="C")
    singleton = satellites.add_satellite("chr1", 2000, 2200, name="D")

    assert satellites.union_by_coordinates(
        "chr1", 90, 210, "chr1", 490, 710
    )
    assert satellites.union_by_coordinates(
        "chr1", 500, 700, "chr1", 1000, 1200
    )

    assert satellites.dsu.find(first) == satellites.dsu.find(second)
    assert satellites.dsu.find(second) == satellites.dsu.find(third)
    assert satellites.dsu.find(singleton) != satellites.dsu.find(first)

    rows = satellites.component_rows()
    linked_rows = rows[:3]
    assert {row["component_id"] for row in linked_rows} == {
        "satellite_component_0001"
    }
    assert {row["component_size"] for row in linked_rows} == {3}
    assert [row["direct_link_count"] for row in linked_rows] == [1, 2, 1]
    assert rows[3]["component_size"] == 1
    assert rows[3]["direct_link_count"] == 0


def test_satellite_dsu_deduplicates_nodes_and_writes_all_members(tmp_path):
    satellites = SatelliteDSU()
    index = satellites.add_satellite("chr2", 10, 50, name=None)
    assert satellites.add_satellite("chr2", 10, 50, name="HSAT") == index

    output_path = tmp_path / "satellite_dsu.tsv"
    satellites.write_tsv(output_path)

    with output_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["satellite_id"] == "satellite_000001"
    assert rows[0]["component_id"] == "satellite_component_0001"
    assert rows[0]["component_size"] == "1"
    assert rows[0]["name"] == "HSAT"

    text_path = tmp_path / "satellite_dsu.txt"
    satellites.write_text(text_path)
    text = text_path.read_text()
    assert "Satellite DSU Results" in text
    assert "Total satellites: 1" in text
    assert "Total components: 1" in text
    assert "[satellite_component_0001] size=1" in text
    assert "chr2:10-50 name=HSAT" in text


def test_satellite_dsu_does_not_union_unmatched_distal_endpoint():
    satellites = SatelliteDSU()
    satellites.add_satellite("chr1", 100, 200)
    satellites.add_satellite("chr1", 500, 600)

    assert not satellites.union_by_coordinates(
        "chr1", 100, 200, "chr1", 800, 900
    )
    assert all(row["component_size"] == 1 for row in satellites.component_rows())


def test_item_rgb_requires_both_dsu_and_matching_ntrprism_values():
    satellites = SatelliteDSU()
    first = satellites.add_satellite("chr1", 100, 200, monomer=171)
    second = satellites.add_satellite("chr1", 300, 400, monomer=171)
    different_prism = satellites.add_satellite("chr1", 500, 600, monomer=5)
    different_dsu = satellites.add_satellite("chr1", 700, 800, monomer=171)
    missing_first = satellites.add_satellite("chr1", 900, 1000, monomer=None)
    missing_second = satellites.add_satellite("chr1", 1100, 1200, monomer=None)
    different_hor = satellites.add_satellite(
        "chr1", 1300, 1400, monomer=171, periodicity=2, is_hor=True
    )
    different_periodicity = satellites.add_satellite(
        "chr1", 1500, 1600, monomer=171, periodicity=[2, 4]
    )
    satellites.union(first, second)
    satellites.union(second, different_prism)
    satellites.union(missing_first, missing_second)
    satellites.union(second, different_hor)
    satellites.union(second, different_periodicity)

    colors = satellites.item_rgb_for_annotations(
        "chr1",
        [100, 300, 500, 700, 900, 1100, 1300, 1500],
        [200, 400, 600, 800, 1000, 1200, 1400, 1600],
        [171, 171, 5, 171, None, None, 171, 171],
        [[2, 3], [2, 3], None, None, None, None, 2, [2, 4]],
        [False, False, False, False, False, False, True, False],
    )

    assert colors[0] == colors[1]
    assert colors[2] != colors[0]
    assert colors[3] != colors[0]
    assert colors[4] != colors[5]
    assert colors[6] != colors[0]
    assert colors[7] != colors[0]
    assert all(
        len(color.split(",")) == 3
        and all(0 <= int(channel) <= 255 for channel in color.split(","))
        for color in colors
    )
