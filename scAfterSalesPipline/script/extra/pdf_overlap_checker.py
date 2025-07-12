import fitz
import math
from PIL import Image
from fitz import Quad, Point
from matplotlib.path import Path
from shapely.geometry import LineString, Polygon
import logging

fitz.TOOLS.set_small_glyph_heights(True)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def rect_dim_len(polygon):
    """计算矩形的长, 宽"""
    if isinstance(polygon, LineString):
        return polygon.length, polygon.length
    elif isinstance(polygon, Polygon):
        x, y = polygon.exterior.coords.xy
        edge_len = (
            Point(x[0], y[0]).distance_to(Point(x[1], y[1])),
            Point(x[1], y[1]).distance_to(Point(x[2], y[2])),
        )
        width = max(edge_len)
        height = min(edge_len)
        return width, height
    else:
        return 0, 0


def polygon_is_rect(polygon):
    """判断多边形是否为矩形"""
    return math.isclose(
        polygon.minimum_rotated_rectangle.area, polygon.area, rel_tol=1e-3
    )


class TextPolygon:
    def __init__(self, page):
        self.upper_chars = ["*", "^", '"', "'"]
        self.lower_chars = ["_", ",", "."]
        self.page = page
        self.blocks = self.page.get_text("rawdict")["blocks"]
        self.mediabox_quad = Quad(
            Point(self.page.mediabox[2], self.page.mediabox[3]),
            Point(self.page.mediabox[0], self.page.mediabox[3]),
            Point(self.page.mediabox[2], self.page.mediabox[1]),
            Point(self.page.mediabox[0], self.page.mediabox[1]),
        )
        self.mediabox_polygon = Polygon(
            [
                self.mediabox_quad.ul,
                self.mediabox_quad.ur,
                self.mediabox_quad.lr,
                self.mediabox_quad.ll,
            ]
        )
        self.text_polygon = self.get_text_polygon()
        self.is_overlap_with_text = False
        self.is_outside = False
        self.is_overlap_with_shape = False
        self.marked_quad = []
        self.shape = self.get_shape()
        self.marked_shape = []

    def get_shape(self):
        """获得图形的坐标"""
        shape = []
        clipped = False
        for i in self.page.get_drawings(extended=True):
            if i.get("scissor"):
                scissor = i.get("scissor")
                clipped = True
            # 忽略线条, 填充颜色为白色/无色/透明的图形是否与文本块相交
            cond1 = (
                (i.get("fill") or (1, 1, 1))
                == (i.get("color") or (1, 1, 1))
                == (1, 1, 1)
            )
            cond2 = (i.get("fill_opacity") or 0) != 1 and (
                i.get("stroke_opacity") or 0
            ) != 1
            if not (cond1 or cond2):
                # 图形具有剪辑层
                if clipped:
                    i["scissor"] = scissor

                i["_shape"] = self.get_shape_polygon(i)
                shape.append(i)
        return shape

    def get_text_polygon(self):
        """获得文本块坐标"""
        text_polygon = []
        for b in self.blocks:
            if b.get("lines"):
                for x in b["lines"]:
                    span = x["spans"][0]
                    if "dingbat" not in span.get("font").lower():
                        text = "".join([c.get("c") for c in span.get("chars")])
                        quad = fitz.recover_quad(line_dir=x["dir"], span=span)
                        if all([c in self.upper_chars for c in set(text)]):
                            hl = (quad.ul + quad.ll) / 2
                            hr = (quad.ur + quad.lr) / 2
                            quad = Quad(quad.ul, quad.ur, hl, hr)
                        elif all([c in self.lower_chars for c in set(text)]):
                            hl = (quad.ul + quad.ll) / 2
                            hr = (quad.ur + quad.lr) / 2
                            quad = Quad(hl, hr, quad.ll, quad.lr)
                        polygon = Polygon([quad.ul, quad.ur, quad.lr, quad.ll])
                        text_polygon.append((text, quad, polygon))
        return text_polygon

    def text_outside(self, **kwargs):
        """文本块是否在页面外"""
        if kwargs.get("ignore_text_outside"):
            return
        for text, quad, polygon in self.text_polygon:
            if not self.mediabox_polygon.contains(polygon):
                self.marked_quad.append(quad)
                self.is_outside = True
        if self.is_outside:
            logger.warning(f"{self.page.parent.name} 发现文本块在页面范围外.")

    def text_overlap(self, **kwargs):
        """文本块之间是否重叠"""
        if kwargs.get("ignore_text_overlap"):
            return
        for i, (text1, quad1, polygon1) in enumerate(self.text_polygon):
            for j, (text2, quad2, polygon2) in enumerate(self.text_polygon[i + 1 :]):
                _f = False
                if polygon1.intersects(polygon2):
                    polygon_inter = polygon1.intersection(polygon2)
                    centroid_inter = polygon_inter.centroid
                    w1, h1 = rect_dim_len(polygon1)
                    w2, h2 = rect_dim_len(polygon2)
                    centroid_d1 = centroid_inter.distance(polygon1.exterior)
                    centroid_d2 = centroid_inter.distance(polygon2.exterior)
                    if all(
                        [
                            polygon_is_rect(p)
                            for p in [polygon1, polygon2, polygon_inter]
                        ]
                    ):
                        wi, hi = rect_dim_len(polygon_inter)
                        # 对于平行的文本块, 相交部分高度超过字体大小的 20% 视为重叠
                        if hi / min(h1, h2) > 0.2:
                            _f, self.is_overlap_with_text = True, True
                    # 重叠区域中心距离文本块边框距离超过字体大小 10% 视为重叠
                    elif centroid_d1 > h1 * 0.1 or centroid_d2 > h2 * 0.1:
                        _f, self.is_overlap_with_text = True, True

                if _f:
                    for i in [quad1, quad2]:
                        if i not in self.marked_quad:
                            self.marked_quad.append(i)

        if self.is_overlap_with_text:
            logger.warning(f"{self.page.parent.name} 发现文本块重叠.")

    def get_shape_polygon(self, shape):
        items = shape.get("items")
        shape_type = items[0][0]
        if shape_type == "c":
            c_codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4] * len(items)
            cubic_bezier = [coord for curve in items for coord in curve[1:]]
            curve_path = Path(cubic_bezier, c_codes)
            coords = [
                coord
                for line in curve_path.to_polygons(closed_only=False)
                for coord in line
            ]
        elif shape_type == "l":
            coords = [coord for line in items for coord in line[1:]]
        elif shape_type == "re":
            quad = Quad(
                items[0][1].top_left,
                items[0][1].top_right,
                items[0][1].bottom_left,
                items[0][1].bottom_right,
            )
            coords = [quad.ul, quad.ur, quad.lr, quad.ll]
        elif shape_type == "qu":
            coords = [items[0][1].ul, items[0][1].ur, items[0][1].lr, items[0][1].ll]

        if shape.get("closePath") or shape_type in ["re", "qu"]:
            try:
                shape_obj = Polygon(coords)
            except Exception:
                shape_obj = LineString(coords)
        else:
            shape_obj = LineString(coords)

        if shape.get("scissor"):
            scissor_quad = Quad(
                shape.get("scissor").top_left,
                shape.get("scissor").top_right,
                shape.get("scissor").bottom_left,
                shape.get("scissor").bottom_right,
            )
            scissor_polygon = Polygon(
                [scissor_quad.ul, scissor_quad.ur, scissor_quad.lr, scissor_quad.ll]
            )
            if shape_obj.intersects(scissor_polygon):
                try:
                    shape_obj = shape_obj.intersection(scissor_polygon)
                except Exception:
                    shape_obj = shape_obj.buffer(1).intersection(scissor_polygon)

        return shape_obj

    def text_shape_overlap(self, ignore=[], **kwargs):
        """文本块与图形是否重叠"""
        if kwargs.get("ignore_shape_overlap"):
            return
        for i, (text, quad, polygon) in enumerate(self.text_polygon):
            # 检查顺序: 顶部 --> 底部
            for shape in self.shape[::-1]:
                _f = False
                if shape.get("items")[0][0] not in ignore:
                    _shape = shape["_shape"]
                    # 跳过点图形, 空图形
                    if _shape.is_empty or (_shape.area == _shape.length == 0):
                        continue
                    if polygon.intersects(_shape):
                        polygon_inter = polygon.intersection(_shape)
                        minrect_inter = polygon_inter.minimum_rotated_rectangle
                        centroid_inter = polygon_inter.centroid
                        wt, ht = rect_dim_len(polygon)
                        wi, hi = rect_dim_len(minrect_inter)
                        if isinstance(_shape, Polygon):
                            cond0 = shape.get("fill") and shape.get("fill_opacity") != 0
                            cond1 = _shape.contains(polygon)
                            cond2 = minrect_inter.equals(polygon)
                            cond3 = min((wt - wi) / ht, (ht - hi) / ht) <= 0.2
                            cond4 = polygon_inter.area / polygon.area > 0.95
                            # 图形"包裹"文本块
                            if cond0 and (cond1 or ((cond2 or cond3) and cond4)):
                                break
                            elif centroid_inter.distance(
                                polygon.exterior
                            ) > ht * 0.1 and not _shape.contains(polygon):
                                _f, self.is_overlap_with_shape = True, True
                        else:
                            if centroid_inter.distance(polygon.exterior) > ht * 0.1:
                                _f, self.is_overlap_with_shape = True, True
                if _f:
                    if _shape not in self.marked_shape:
                        self.marked_shape.append(shape.get("items"))
                    if quad not in self.marked_quad:
                        self.marked_quad.append(quad)

        if self.is_overlap_with_shape:
            logger.warning(f"{self.page.parent.name} 发现文本块与图形重叠.")

    def draw_page(self):
        """重叠部分高亮显示"""
        for i in self.marked_quad:
            self.page.draw_quad(
                i,
                color=(1, 0, 0),
                width=1,
                overlay=True,
                fill=(1, 0, 0),
                fill_opacity=0.1,
            )
        for i in self.marked_shape:
            for item in i:
                if item[0] == "qu":
                    self.page.draw_quad(
                        item[1],
                        color=(1, 0, 0),
                        width=1,
                        overlay=True,
                        fill=(1, 0, 0),
                        fill_opacity=0,
                    )
                elif item[0] == "l":
                    self.page.draw_line(
                        *item[1:], color=(1, 0, 0), width=1, overlay=True
                    )
                elif item[0] == "c":
                    self.page.draw_bezier(
                        *item[1:], color=(1, 0, 0), width=1, overlay=True
                    )
                elif item[0] == "re":
                    self.page.draw_rect(
                        item[1],
                        color=(1, 0, 0),
                        width=1,
                        overlay=True,
                        fill=(1, 0, 0),
                        fill_opacity=0,
                    )

        pix = self.page.get_pixmap(dpi=100)
        mode = "RGBA" if pix.alpha else "RGB"
        img = Image.frombytes(mode, [pix.width, pix.height], pix.samples)
        return img

    def check(self, **kwargs):
        ignore_shape_type = [
            (kwargs.get("ignore_line"), "l"),
            (kwargs.get("ignore_curve"), "c"),
            (kwargs.get("ignore_rect"), "re"),
            (kwargs.get("ignore_quad"), "qu"),
        ]
        ignore = [j for i, j in ignore_shape_type if i]
        logging.disable(logging.WARNING)
        self.text_outside(**kwargs)
        self.text_overlap(**kwargs)
        self.text_shape_overlap(ignore=ignore, **kwargs)
        logging.disable(logging.NOTSET)
        check_result = [
            (self.is_outside, "文本块在页面范围外"),
            (self.is_overlap_with_text, "文本块重叠"),
            (self.is_overlap_with_shape, "文本块与图形重叠"),
        ]
        message = ", ".join([j for i, j in check_result if i])
        if message:
            logger.warning(f"{self.page.parent.name} 发现{message}.")


if __name__ == "__main__":
    import click
    import sys
    import tempfile
    import pathlib

    @click.command("PDF overlap checker")
    @click.argument("files", type=click.Path(exists=True), nargs=-1)
    @click.option("--ignore_outside", is_flag=True, help="忽略页面外的文本块")
    @click.option("--ignore_text_overlap", is_flag=True, help="忽略重叠的文本块")
    @click.option("--ignore_shape_overlap", is_flag=True, help="忽略与图形重叠的文本块")
    @click.option("--ignore_line", is_flag=True, help="忽略线段的重叠检查")
    @click.option("--ignore_curve", is_flag=True, help="忽略曲线的重叠检查")
    @click.option("--ignore_rect", is_flag=True, help="忽略矩形的重叠检查")
    @click.option("--ignore_quad", is_flag=True, help="忽略多边形的重叠检查")
    @click.option(
        "--draw",
        default=None,
        type=click.Path(exists=False, file_okay=False),
        help="异常图片输出目录",
    )
    def main(files, **kwargs):
        # print(kwargs)
        if not sys.stdin.isatty():
            files = [line.strip() for line in sys.stdin]
            sys.argv += files

        draw = kwargs.get("draw")
        if draw:
            warn_dir = draw
            pathlib.Path(warn_dir).mkdir(parents=True, exist_ok=True)
            # warn_dir = tempfile.mkdtemp(prefix='pdf-overlap.warning.',
            #                             dir=pathlib.Path.cwd())
            logger.info(f"图片目录: {warn_dir}")

        for file in files:
            with fitz.open(file) as doc:
                page = doc.load_page(0)
                tp = TextPolygon(page)
                tp.check(**kwargs)

                if draw:
                    isOverlapped = any(
                        [
                            tp.is_outside,
                            tp.is_overlap_with_text,
                            tp.is_overlap_with_shape,
                        ]
                    )
                    if isOverlapped:
                        prefix = pathlib.Path(file).stem
                        suffix = next(tempfile._get_candidate_names())
                        draw_name = f"{prefix}.{suffix}.png"
                        img = tp.draw_page()
                        img.save(str(pathlib.Path(warn_dir, draw_name)))

    main()
