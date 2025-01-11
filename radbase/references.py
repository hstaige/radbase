from dataclasses import dataclass
from typing import Optional


@dataclass
class Reference:
    title: str
    author: list[str]
    year: str
    month: Optional[str] = None
    journal: Optional[str] = None
    pages: Optional[int] = None
    volume: Optional[int] = None
    url: Optional[str] = None
    doi: Optional[str] = None
